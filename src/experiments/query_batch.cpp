#include "experiments/experiment.h"
#include "common/memory_guard.h"
#include "index/index_io.h"
#include "index/index_space.h"
#include "query/index_query_utils.h"
#include "query/mat_recovery.h"
#include "query/online_query.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <sys/resource.h>
#include <unistd.h>

namespace tacc {
namespace {
namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
constexpr std::size_t query_count = 500, empty_quota = 450, nonempty_quota = 50;
constexpr std::uint64_t seed = 20261023;
using Key = std::tuple<Timestamp,Timestamp,Count,Count,AttributeMask,AttributeMask>;
using Fields = std::map<std::string,std::string>;
void require(bool ok, const std::string& why) { if (!ok) throw std::runtime_error(why); }
double elapsed_us(Clock::time_point start) {
    return std::chrono::duration<double,std::micro>(Clock::now()-start).count();
}
AttributeMask mask(const std::vector<AttributeID>& values) {
    AttributeMask m=0;
    for (auto v:values) { require(v<64,"attribute ID exceeds 63"); m|=AttributeMask{1}<<v; }
    return m;
}
std::vector<AttributeID> attributes(AttributeMask m) {
    std::vector<AttributeID> out;
    for (AttributeID a=0;a<64;++a) if (m&(AttributeMask{1}<<a)) out.push_back(a);
    return out;
}
Key key(const Query& q) { return {q.window.start,q.window.end,q.alpha,q.beta,mask(q.u_attributes),mask(q.v_attributes)}; }
Fields fields(const std::string& line) {
    Fields out; std::istringstream in(line); std::string token;
    while (in>>token) { const auto n=token.find('='); if (n!=std::string::npos) out[token.substr(0,n)]=token.substr(n+1); }
    return out;
}
std::string first_line(const fs::path& file) {
    std::ifstream in(file); std::string s;
    require(bool(std::getline(in,s)),"missing/empty file: "+file.string()); return s;
}
std::uint64_t number(const std::string& s) {
    require(!s.empty() && std::all_of(s.begin(),s.end(),[](char c){return c>='0'&&c<='9';}),"invalid unsigned integer");
    return std::stoull(s);
}
std::uint64_t query_hash(const std::vector<Query>& queries) {
    std::uint64_t h=14695981039346656037ULL;
    for (const auto& q:queries) for (std::uint64_t x:{std::uint64_t(q.window.start),std::uint64_t(q.window.end),std::uint64_t(q.alpha),std::uint64_t(q.beta),mask(q.u_attributes),mask(q.v_attributes)})
        for (int i=0;i<8;++i) { h^=x&255;h*=1099511628211ULL;x>>=8; }
    return h;
}
std::uint64_t file_hash(const fs::path& path) {
    std::ifstream in(path,std::ios::binary); require(bool(in),"cannot hash file: "+path.string());
    std::array<char,65536> buffer{}; std::uint64_t h=14695981039346656037ULL;
    while (in) { in.read(buffer.data(),buffer.size()); for (std::streamsize i=0;i<in.gcount();++i) {h^=static_cast<unsigned char>(buffer[i]);h*=1099511628211ULL;} }
    require(in.eof(),"cannot hash complete file"); return h;
}
std::uint64_t executable_hash() { return file_hash("/proc/self/exe"); }
std::string host_name() {
    std::array<char,256> name{};require(gethostname(name.data(),name.size()-1)==0,"cannot identify host");return name.data();
}
TemporalBipartiteGraph load_graph(const fs::path& data) {
    TemporalBipartiteGraph graph;
    require(graph.load_from_files((data/"graph_ulabel.txt").string(),(data/"graph_vlabel.txt").string(),(data/"graph_timestamp.txt").string()),"graph load failed");
    return graph;
}
AttributeMask universe(const TemporalBipartiteGraph& g,bool u) {
    AttributeMask out=0;
    for (VertexID v=0;v<(u?g.num_u_vertices():g.num_v_vertices());++v) {
        const auto a=u?g.u_attribute(v):g.v_attribute(v); require(a<64,"unsupported attributes");out|=AttributeMask{1}<<a;
    }
    return out;
}
bool empty(const QueryResult& r) {
    require(r.u_vertices.empty()==r.v_vertices.empty(),"invalid one-sided core result");
    return r.u_vertices.empty();
}
void sort_result(QueryResult& r) {
    std::sort(r.u_vertices.begin(),r.u_vertices.end());std::sort(r.v_vertices.begin(),r.v_vertices.end());
}
void metadata(std::ostream& out,const TemporalBipartiteGraph& g,std::uint64_t hash) {
    out<<"META format=UV4501 graph_fingerprint="<<graph_fingerprint(g)<<" queries=500 target_empty=450 target_nonempty=50"
       <<" result_selection=class_quota query_hash="<<hash<<" seed="<<seed<<" rounds=1 warmup=0\n";
}
void check_metadata(const std::string& line,const TemporalBipartiteGraph& graph,std::uint64_t hash=0) {
    const auto f=fields(line);
    require(f.at("graph_fingerprint")==std::to_string(graph_fingerprint(graph)) && f.at("queries")=="500" &&
        f.at("target_empty")=="450" && f.at("target_nonempty")=="50" && f.at("result_selection")=="class_quota","workload metadata mismatch");
    if (hash) require(f.at("query_hash")==std::to_string(hash),"query condition hash mismatch");
}
std::vector<Query> read_queries(const fs::path& path,const TemporalBipartiteGraph& graph) {
    std::ifstream in(path); std::string line;
    require(bool(std::getline(in,line)),"cannot read conditions");check_metadata(line,graph);
    const auto meta=fields(line);
    require(bool(std::getline(in,line)) && line=="id\tts\tte\talpha\tbeta\tu_mask\tv_mask","invalid condition header");
    std::vector<Query> out;std::set<Key> seen;
    const auto um=universe(graph,true),vm=universe(graph,false);
    while (std::getline(in,line)) {
        std::istringstream row(line);std::uint64_t id,s,e,a,b,u,v;std::string extra;
        require(bool(row>>id>>s>>e>>a>>b>>u>>v) && !(row>>extra) && id==out.size()+1 && id<=query_count,"invalid condition row");
        require(s<=e && e<=std::numeric_limits<Timestamp>::max() && a>0 && a<=std::min<std::size_t>(10,graph.num_v_vertices()) &&
            b>0 && b<=std::min<std::size_t>(10,graph.num_u_vertices()) && u && v && !(u&~um) && !(v&~vm),"condition out of domain");
        require(std::binary_search(graph.time_points().begin(),graph.time_points().end(),static_cast<Timestamp>(s))&&
            std::binary_search(graph.time_points().begin(),graph.time_points().end(),static_cast<Timestamp>(e)),"condition endpoint is not an observed time point");
        Query q{{static_cast<Timestamp>(s),static_cast<Timestamp>(e)},attributes(u),attributes(v),static_cast<Count>(a),static_cast<Count>(b)};
        require(seen.insert(key(q)).second,"duplicate condition");out.push_back(std::move(q));
    }
    require(out.size()==query_count,"need exactly 500 conditions");
    if (meta.count("query_hash")) require(meta.at("query_hash")==std::to_string(query_hash(out)),"condition file corrupted");
    return out;
}
struct Discovery { std::size_t draws{},calls{},rejected{}; double seconds{}; };
std::vector<Query> generate(const TemporalBipartiteGraph& g,Discovery& audit) {
    const auto started=Clock::now();
    const auto uu=attributes(universe(g,true)),vv=attributes(universe(g,false));
    require(!uu.empty()&&!vv.empty()&&uu.size()<64&&vv.size()<64,"invalid attribute universe");
    const auto subset=[](const auto& values,AttributeMask m) { std::vector<AttributeID> out;for (std::size_t i=0;i<values.size();++i) if(m&(AttributeMask{1}<<i))out.push_back(values[i]);return out; };
    std::mt19937_64 random(seed);
    const auto& times=g.time_points();std::uniform_int_distribution<std::size_t> pick_start(0,times.size()-1);
    std::uniform_int_distribution<Count> pick_a(1,std::min<std::size_t>(10,g.num_v_vertices())),pick_b(1,std::min<std::size_t>(10,g.num_u_vertices()));
    std::uniform_int_distribution<AttributeMask> pick_u(1,(AttributeMask{1}<<uu.size())-1),pick_v(1,(AttributeMask{1}<<vv.size())-1);
    std::vector<Query> out;std::set<Key> seen;std::size_t ne=0,nn=0;
    for (;out.size()<query_count && audit.draws<1000000;++audit.draws) {
        const auto s=pick_start(random),e=std::uniform_int_distribution<std::size_t>(s,times.size()-1)(random);
        const auto a=pick_a(random),b=pick_b(random);
        const auto u=pick_u(random),v=pick_v(random);
        Query q{{times[s],times[e]},subset(uu,u),subset(vv,v),a,b};
        if (!seen.insert(key(q)).second) continue;
        const auto r=online_query(g,q);++audit.calls;
        const bool is_empty=empty(r);auto& n=is_empty?ne:nn;
        if (n<(is_empty?empty_quota:nonempty_quota)) {++n;out.push_back(std::move(q));} else ++audit.rejected;
        if (audit.calls%500==0) std::cout<<"PROGRESS phase=quota_discovery calls="<<audit.calls<<" empty="<<ne<<" nonempty="<<nn<<std::endl;
    }
    require(ne==empty_quota && nn==nonempty_quota,"cannot find 450 distinct empty + 50 distinct nonempty conditions in fixed domain; no quota/threshold relaxation performed");
    std::shuffle(out.begin(),out.end(),random);audit.seconds=elapsed_us(started)/1e6;return out;
}
void write_ids(std::ostream& out,const std::vector<VertexID>& ids) {
    if (ids.empty()) {out<<'-';return;} for (std::size_t i=0;i<ids.size();++i){if(i)out<<',';out<<ids[i];}
}
std::vector<VertexID> read_ids(const std::string& s) {
    if(s=="-")return {};
    std::vector<VertexID> out;std::istringstream in(s);std::string v;
    while(std::getline(in,v,',')){const auto x=number(v);require(x>0&&x<=std::numeric_limits<VertexID>::max(),"bad answer ID");out.push_back(static_cast<VertexID>(x));}
    require(!out.empty() && std::is_sorted(out.begin(),out.end()) && std::adjacent_find(out.begin(),out.end())==out.end(),"invalid answer order");return out;
}
struct Total { std::size_t count{};double us{}; };
using Totals=std::array<Total,3>;
void add(Totals& totals,bool is_empty,double us) {for (auto i:{0u,is_empty?1u:2u}){++totals[i].count;totals[i].us+=us;}}
void summaries(std::ostream& out,const char* engine,const Totals& totals,const Totals* baseline=nullptr,double prep_us=0) {
    const std::array<const char*,3> names{{"all","empty","nonempty"}};
    for(std::size_t i=0;i<3;++i){const auto& t=totals[i];out<<"SUMMARY engine="<<engine<<" group="<<names[i]<<" queries="<<t.count<<" total_us="<<t.us<<" avg_us="<<(t.count?t.us/t.count:0);
        if(baseline)out<<" speedup_vs_online="<<(t.us?(*baseline)[i].us/t.us:0);
        if(i==0&&baseline)out<<" amortized_avg_us="<<(t.us+prep_us)/query_count<<" speedup_including_aux="<<((*baseline)[i].us/(t.us+prep_us));
        out<<'\n';}
}
void close_checked(std::ofstream& out) { out.close();require(bool(out),"output write/close failed"); }

bool batch_query_profile_enabled() {
    const char* value=std::getenv("TAC_QUERY_PROFILE");
    if(value==nullptr || std::string(value)=="0")return false;
    if(std::string(value)=="1")return true;
    throw std::invalid_argument("TAC_QUERY_PROFILE must be 0 or 1");
}

void online_stage(const fs::path& data,const fs::path& conditions,const fs::path& output) {
    auto graph=load_graph(data);require(fs::create_directory(output),"Online output directory must be new");
    const bool existing=fs::exists(conditions);Discovery audit;
    auto queries=existing?read_queries(conditions,graph):generate(graph,audit);
    const auto hash=query_hash(queries),binary_hash=executable_hash();
    std::ofstream report(output/"online.log"),samples(output/"online.samples.csv"),answers(output/"answers.tsv");
    require(bool(report)&&bool(samples)&&bool(answers),"cannot create Online outputs");
    report<<std::fixed<<std::setprecision(6);samples<<std::fixed<<std::setprecision(6);
    metadata(report,graph,hash);metadata(answers,graph,hash);metadata(samples,graph,hash);
    report<<"DISCOVERY source="<<(existing?"loaded_conditions":"generated")<<" draws="<<audit.draws<<" online_calls="<<audit.calls<<" discovery_s="<<audit.seconds<<" rejected_quota="<<audit.rejected<<" selection_by_timing=0 discovery_excluded_from_query_time=1\n";
    samples<<"id,empty,online_us,u_vertices,v_vertices\n";answers<<"id\tu_ids\tv_ids\n";
    Totals totals{};
    for(std::size_t i=0;i<queries.size();++i){
        const auto start=Clock::now();auto result=online_query(graph,queries[i]);const double duration=elapsed_us(start);
        const bool is_empty=empty(result);add(totals,is_empty,duration);sort_result(result);
        samples<<i+1<<','<<is_empty<<','<<duration<<','<<result.u_vertices.size()<<','<<result.v_vertices.size()<<'\n';
        answers<<i+1<<'\t';write_ids(answers,result.u_vertices);answers<<'\t';write_ids(answers,result.v_vertices);answers<<'\n';
        if((i+1)%100==0)std::cout<<"PROGRESS phase=quota_online completed="<<i+1<<std::endl;
    }
    require(totals[1].count==empty_quota&&totals[2].count==nonempty_quota,"fresh Online does not match required quota");
    summaries(report,"online",totals);summaries(std::cout,"online",totals);
    report<<"RESULT exp=uv_batch engine=online queries=500 empty=450 nonempty=50 rounds=1 timing=fresh status=ok\n";
    close_checked(report);close_checked(samples);close_checked(answers);
    if(!existing){
        if(!conditions.parent_path().empty())fs::create_directories(conditions.parent_path());
        const auto pending=conditions.string()+".pending."+std::to_string(getpid());
        std::ofstream out(pending);require(bool(out),"cannot write conditions");metadata(out,graph,hash);out<<"id\tts\tte\talpha\tbeta\tu_mask\tv_mask\n";
        for(std::size_t i=0;i<queries.size();++i){const auto& q=queries[i];out<<i+1<<'\t'<<q.window.start<<'\t'<<q.window.end<<'\t'<<q.alpha<<'\t'<<q.beta<<'\t'<<mask(q.u_attributes)<<'\t'<<mask(q.v_attributes)<<'\n';}
        close_checked(out);require(!fs::exists(conditions),"conditions appeared concurrently; use the workflow lock");fs::rename(pending,conditions);
    }
    std::ofstream complete(output/"online.complete");require(bool(complete),"cannot write completion marker");
    complete<<"COMPLETE status=ok graph_fingerprint="<<graph_fingerprint(graph)<<" query_hash="<<hash<<" binary_hash="<<binary_hash
        <<" host="<<host_name()<<" answers_hash="<<file_hash(output/"answers.tsv")<<" samples_hash="<<file_hash(output/"online.samples.csv")
        <<" queries=500 empty=450 nonempty=50 timing=fresh\n";close_checked(complete);
    std::cout<<"RESULT exp=uv_batch engine=online queries=500 empty=450 nonempty=50 rounds=1 timing=fresh status=ok\n";
}

template<class Index> QueryResult ask(const Index& index,const TemporalBipartiteGraph& graph,const Query& q,QueryProfile* p) {
    if constexpr(std::is_same_v<Index,MATIndex>)return query_mat_block_index(index,graph,q,p);
    else if constexpr(std::is_same_v<Index,UATIndex>)return query_uat_block_index(index,graph,q,p);
    else return query_membership_index(index,graph,q,p);
}
template<class Index> void query_stage(const std::string& kind,const fs::path& data,const fs::path& binary,
                                      const fs::path& conditions,const fs::path& online,const fs::path& output) {
    const bool collect_profile=batch_query_profile_enabled();
    const bool uses_grouped=index_uses_grouped_recovery(kind);
    const auto bytes=fs::file_size(binary);
    const auto mtime=fs::last_write_time(binary);
    auto graph=load_graph(data);const auto queries=read_queries(conditions,graph);const auto hash=query_hash(queries);
    const auto complete=fields(first_line(online/"online.complete"));
    require(complete.at("status")=="ok"&&complete.at("graph_fingerprint")==std::to_string(graph_fingerprint(graph))&&complete.at("query_hash")==std::to_string(hash)&&
        complete.at("binary_hash")==std::to_string(executable_hash())&&complete.at("queries")=="500"&&complete.at("empty")=="450"&&complete.at("nonempty")=="50"&&complete.at("timing")=="fresh","Online run is incomplete or belongs to different data/conditions/executable");
    require(complete.at("host")==host_name()&&complete.at("answers_hash")==std::to_string(file_hash(online/"answers.tsv"))&&
        complete.at("samples_hash")==std::to_string(file_hash(online/"online.samples.csv")),"Online host or answer/timing file integrity mismatch");
    IndexScope scope;std::string error;require(read_index_scope(binary.string(),scope,&error)&&scope==IndexScope::UV,"expected UV binary: "+error);
    Index index;auto started=Clock::now();require(load_index_binary(binary.string(),graph,index,&error),"index load failed: "+error);const auto load_us=elapsed_us(started);
    require(fs::file_size(binary)==bytes&&fs::last_write_time(binary)==mtime,"index changed during loading");
    require(index.stats().method==BuildMethod::Opt,"expected opt index");
    const auto before_hash=formal_index_hash(index);
    require(fs::create_directory(output),"index output directory must be new");
    std::ofstream report(output/(kind+".log")),samples(output/(kind+".samples.csv"));require(bool(report)&&bool(samples),"cannot create query outputs");
    report<<std::fixed<<std::setprecision(6);samples<<std::fixed<<std::setprecision(6);metadata(report,graph,hash);
    report<<"MEASUREMENT query_profile="<<collect_profile<<" timing=outer_wall_clock"
          <<" diagnostic_columns="<<(collect_profile?"measured":"NA")<<'\n';
    double prep_us=0;
    if(uses_grouped){started=Clock::now();const auto& aux=graph.grouped_v_time_adjacency();prep_us=elapsed_us(started);
            report<<"RESULT exp="<<kind<<"_query_aux prepare_us="<<prep_us<<" aux_bytes="<<aux.bytes()<<" scratch_peak_bytes="<<aux.scratch_peak_bytes<<" groups="<<aux.times.size()<<" events="<<aux.neighbors.size()<<" status=ok\n";}
    const auto space=formal_index_space(index);
    write_formal_index_space(report,kind.c_str(),space,index.stats().approx_bytes);
    write_query_index_space(report,kind.c_str(),space,graph,uses_grouped);
    report<<"INDEX kind="<<kind<<" scope=uv query_route=uv source=loaded built=0 hash="<<before_hash<<" load_us="<<load_us<<" aux_prepare_us="<<prep_us
        <<" mat_recovery="<<(kind=="mat"?mat_recovery_name(configured_mat_recovery()):"not_applicable")
        <<" uat_recovery="<<(kind=="uat"?uat_recovery_name(configured_uat_recovery()):"not_applicable")
        <<" seed_recovery="<<seed_recovery_name(configured_seed_recovery())<<" rounds=1 warmup=0 path="<<binary.string()<<'\n';
    samples<<"id,empty,index_us,online_us,same,window_us,peeling_us,clear_us,window_edges,grouped_events,retained_bytes\n";
    std::ifstream answers(online/"answers.tsv"),baseline(online/"online.samples.csv");std::string line;
    require(bool(std::getline(answers,line)),"missing answers metadata");check_metadata(line,graph,hash);
    require(bool(std::getline(answers,line))&&line=="id\tu_ids\tv_ids","invalid answers header");
    require(bool(std::getline(baseline,line)),"missing baseline metadata");check_metadata(line,graph,hash);
    require(bool(std::getline(baseline,line))&&line=="id,empty,online_us,u_vertices,v_vertices","invalid baseline header");
    Totals totals{},online_totals{};QueryProfile profile;std::size_t mismatches=0;
    for(std::size_t i=0;i<queries.size();++i){
        require(bool(std::getline(answers,line)),"truncated Online answers");std::istringstream row(line);std::string id,u,v,extra;
        require(bool(row>>id>>u>>v)&&!(row>>extra)&&number(id)==i+1,"invalid answer row");QueryResult reference{read_ids(u),read_ids(v)};
        require(bool(std::getline(baseline,line)),"truncated Online timing");std::replace(line.begin(),line.end(),',',' ');std::istringstream br(line);
        std::size_t bid,nu,nv;int is_empty;double base;
        require(bool(br>>bid>>is_empty>>base>>nu>>nv)&&!(br>>extra)&&bid==i+1&&nu==reference.u_vertices.size()&&nv==reference.v_vertices.size()&&is_empty==int(empty(reference))&&std::isfinite(base)&&base>=0,"Online timing/reference mismatch");
        QueryProfile current;started=Clock::now();auto actual=ask(index,graph,queries[i],collect_profile?&current:nullptr);const auto duration=elapsed_us(started);sort_result(actual);
        const bool same=actual.u_vertices==reference.u_vertices&&actual.v_vertices==reference.v_vertices;mismatches+=!same;
        add(totals,is_empty,duration);add(online_totals,is_empty,base);
        if(collect_profile)profile+=current;
        samples<<i+1<<','<<is_empty<<','<<duration<<','<<base<<','<<same;
        if(collect_profile)samples<<','<<current.window_recovery_us<<','<<current.peeling_us<<','<<current.workspace_clear_us<<','<<current.window_edges<<','<<current.grouped_window_events<<','<<current.workspace_retained_bytes;
        else samples<<",NA,NA,NA,NA,NA,NA";
        samples<<'\n';
        if((i+1)%100==0)std::cout<<"PROGRESS phase=quota_query engine="<<kind<<" completed="<<i+1<<std::endl;
    }
    require(!std::getline(answers,line)&&!std::getline(baseline,line),"unexpected trailing answer/timing rows");
    require(totals[1].count==empty_quota&&totals[2].count==nonempty_quota,"incorrect result-class quota");
    require(formal_index_hash(index)==before_hash,"query mutated index");
    summaries(report,kind.c_str(),totals,&online_totals,prep_us);summaries(std::cout,kind.c_str(),totals,&online_totals,prep_us);
    if(collect_profile)report<<"PROFILE enabled=1 grouped_recoveries="<<profile.grouped_recoveries<<" grouped_events="<<profile.grouped_window_events<<" budget_fallbacks="<<profile.grouped_budget_fallbacks
        <<" workspace_fixed_bytes="<<profile.workspace_fixed_bytes<<" retained_bytes="<<profile.workspace_retained_bytes;
    else report<<"PROFILE enabled=0 counters=not_collected";
    struct rusage usage {};
    if(getrusage(RUSAGE_SELF,&usage)==0) report<<" peak_rss_status=ok hwm_kb="<<usage.ru_maxrss
        <<" peak_rss_mb="<<space_mb(std::uint64_t(usage.ru_maxrss)*1024)
        <<" peak_rss_G="<<space_G(std::uint64_t(usage.ru_maxrss)*1024)
        <<" peak_scope=process_through_query\n";
    else report<<" peak_rss_status=unavailable hwm_kb=NA\n";
    report<<"RESULT exp=uv_batch engine="<<kind<<" source=loaded built=0 queries=500 empty=450 nonempty=50 rounds=1 mismatches="<<mismatches<<" query_hash="<<hash<<" status="<<(mismatches?"failed":"ok")<<'\n';
    close_checked(report);close_checked(samples);require(mismatches==0,"index differs from fresh Online answers");
    std::cout<<"RESULT exp=uv_batch engine="<<kind<<" source=loaded built=0 queries=500 empty=450 nonempty=50 rounds=1 mismatches=0 status=ok\n";
}
}

int run_quota_experiment(int argc,char** argv) try {
    const auto guard=configure_memory_guard();write_memory_guard(std::cerr,guard);require(memory_guard_ready(guard),"memory guard failed");
    const std::string command=argc>1?argv[1]:"";
    require((command=="quota-online"&&argc==5)||(command=="quota-query"&&argc==8),
        "Usage: quota-online DATA CONDITIONS_TSV NEW_ONLINE_DIR | quota-query KIND DATA INDEX CONDITIONS_TSV ONLINE_DIR NEW_QUERY_DIR");
    if(command=="quota-online")online_stage(argv[2],argv[3],argv[4]);
    else {
        const std::string kind=argv[2];
        require(configured_query_route()==QueryRoute::Auto||configured_query_route()==QueryRoute::UV,"quota batch requires UV/auto query route");
        if(kind=="mat")query_stage<MATIndex>(kind,argv[3],argv[4],argv[5],argv[6],argv[7]);
        else if(kind=="uat")query_stage<UATIndex>(kind,argv[3],argv[4],argv[5],argv[6],argv[7]);
        else if(kind=="vat")query_stage<VATIndex>(kind,argv[3],argv[4],argv[5],argv[6],argv[7]);
        else throw std::runtime_error("unknown index kind");
    }
    return 0;
}catch(const std::exception& e){std::cerr<<"RESULT exp=uv_batch status=failed reason="<<e.what()<<'\n';return 2;}
}
