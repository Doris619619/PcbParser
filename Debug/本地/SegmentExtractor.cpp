//#include <iostream>
//#include <fstream>
//#include <string>
//#include <vector>
//#include <cmath>
//#include <map>
//#include <unordered_map>
//#include <algorithm>
//#include <cctype>
//#include <stdexcept>
//#include "blocked_area_analyzer.h"  
//// ====================== 基础结构 ======================
//struct Segment {
//    double startX = 0, startY = 0;
//    double endX = 0, endY = 0;
//    double width = 0;
//    std::string layer;
//    int net = 0;
//};
//struct SegmentSpan { size_t begin = 0, end = 0; };
//
//static std::string ReadWholeFile(const std::string& path) {
//    std::ifstream fin(path, std::ios::binary);
//    if (!fin) throw std::runtime_error("无法打开: " + path);
//    fin.seekg(0, std::ios::end);
//    std::string s(fin.tellg(), '\0');
//    fin.seekg(0, std::ios::beg);
//    fin.read(&s[0], s.size());
//    return s;
//}
//static void WriteWholeFile(const std::string& path, const std::string& c) {
//    std::ofstream f(path, std::ios::binary); if (!f) throw std::runtime_error("写出失败"); f << c;
//}
//static const char* skip_ws(const char* p) { while (*p && isspace((unsigned char)*p))++p; return p; }
//static const char* parse_double(const char* p, double& v) { p = skip_ws(p); char* e = nullptr; v = strtod(p, &e); return e == p ? nullptr : e; }
//static const char* parse_int(const char* p, int& v) { p = skip_ws(p); char* e = nullptr; v = strtol(p, &e, 10); return e == p ? nullptr : e; }
//static size_t findNextSexpHead(const std::string& t, size_t from, const std::string& kw) {
//    for (size_t i = from; i < t.size(); ++i)if (t[i] == '(') {
//        size_t j = i + 1; while (j < t.size() && isspace((unsigned char)t[j]))++j;
//        if (t.compare(j, kw.size(), kw) == 0)return i;
//    }return std::string::npos;
//}
//static const char* find_kw(const std::string& s, const char* kw) {
//    const char* t = s.c_str(); for (size_t i = 0; i < s.size(); ++i) {
//        if (t[i] != '(')continue; size_t j = i + 1; while (isspace((unsigned char)t[j]))++j;
//        const char* pj = t + j, * pk = kw; while (*pk && *pj == *pk) { ++pj; ++pk; }if (!*pk)return pj;
//    }return nullptr;
//}
//static bool parse_layer_name(const char* p, std::string& o) {
//    p = skip_ws(p); if (*p == '"') { ++p; const char* q = p; while (*q && *q != '"')++q; o.assign(p, q); }
//    else { const char* q = p; while (*q && !isspace((unsigned char)*q) && *q != ')')++q; o.assign(p, q); }return true;
//}
//static inline void normalizeLayer(std::string& L) { if (L == "Bottom")L = "B.Cu"; else if (L == "Top")L = "F.Cu"; }
//static bool ParseSegmentBlock_NoRegex(const std::string& blk, Segment& s) {
//    const char* p;
//    if ((p = find_kw(blk, "start"))) { p = parse_double(p, s.startX); p = parse_double(p, s.startY); }
//    else return false;
//    if ((p = find_kw(blk, "end"))) { p = parse_double(p, s.endX); p = parse_double(p, s.endY); }
//    else return false;
//    if ((p = find_kw(blk, "width"))) { p = parse_double(p, s.width); }
//    else return false;
//    if ((p = find_kw(blk, "layer"))) { parse_layer_name(p, s.layer); normalizeLayer(s.layer); }
//    else return false;
//    if ((p = find_kw(blk, "net"))) { parse_int(p, s.net); }
//    else return false; return true;
//}
//static void LoadSegmentsWithSpans(const std::string& pcb, std::vector<Segment>& outS, std::vector<SegmentSpan>& outSp, std::string* raw = nullptr) {
//    std::string t = ReadWholeFile(pcb); if (raw)*raw = t; size_t i = 0; while (true) {
//        size_t s = findNextSexpHead(t, i, "segment"); if (s == std::string::npos)break;
//        int dep = 0; size_t j = s; bool st = false; for (; j < t.size(); ++j) { if (t[j] == '(') { ++dep; st = true; } else if (t[j] == ')' && st && --dep == 0) { ++j; break; } }
//        std::string blk = t.substr(s, j - s); Segment seg; if (ParseSegmentBlock_NoRegex(blk, seg)) { outS.push_back(seg); outSp.push_back({ s,j }); }i = j;
//    }
//}
//
//// ====================== 几何计算 ======================
//static inline double sq(double x) { return x * x; }
//static inline double dist2(double x1, double y1, double x2, double y2) { return sq(x1 - x2) + sq(y1 - y2); }
//static double pointSegmentDist2(double px, double py, double ax, double ay, double bx, double by, bool& on) {
//    double vx = bx - ax, vy = by - ay, wx = px - ax, wy = py - ay, vv = vx * vx + vy * vy, t = (vv > 0) ? (wx * vx + wy * vy) / vv : 0; on = true;
//    if (t < 0) { on = false; return dist2(px, py, ax, ay); }if (t > 1) { on = false; return dist2(px, py, bx, by); }double cx = ax + t * vx, cy = ay + t * vy; return dist2(px, py, cx, cy);
//}
//struct MatchParams { double coordTol = 0.05; bool requireLayerEqual = true; };
//static bool endpointTouchesOrig(const Segment& s, double px, double py, const Segment& t, const MatchParams& mp) {
//    if (mp.requireLayerEqual && s.layer != t.layer)return false; bool on; double d2 = pointSegmentDist2(px, py, t.startX, t.startY, t.endX, t.endY, on);
//    return on && d2 <= mp.coordTol * mp.coordTol;
//}
//
//// ====================== 匹配与聚类 ======================
//static void collectHits(const Segment& s, const std::vector<Segment>& orig, const MatchParams& mp, std::vector<std::pair<int, double>>& out) {
//    out.clear(); for (auto& t : orig) {
//        if (mp.requireLayerEqual && s.layer != t.layer)continue;
//        if (endpointTouchesOrig(s, s.startX, s.startY, t, mp))out.emplace_back(t.net, t.width);
//        if (endpointTouchesOrig(s, s.endX, s.endY, t, mp))out.emplace_back(t.net, t.width);
//    }
//}
//static void pickOrigNetAndWidthForSmallSegment(const Segment& s, const std::vector<Segment>& orig, const MatchParams& mp,
//    int& net, double& w, std::vector<std::pair<int, double>>* outHits = nullptr) {
//    std::vector<std::pair<int, double>> hits; collectHits(s, orig, mp, hits); if (outHits)*outHits = hits;
//    if (hits.empty()) { net = -1; w = s.width; return; }
//    std::map<int, int>vNet; for (auto& h : hits)vNet[h.first]++; net = -1; int bv = -1; for (auto& kv : vNet)
//        if (kv.second > bv || (kv.second == bv && kv.first < net)) { bv = kv.second; net = kv.first; }
//    std::map<double, int>vW; for (auto& h : hits)if (h.first == net)vW[h.second]++; int wv = -1; w = 0; for (auto& kv : vW)
//        if (kv.second > wv || (kv.second == wv && kv.first < w)) { wv = kv.second; w = kv.first; }
//}
//
//struct ClusterNode { double x = 0, y = 0; std::vector<int>segs; };
//static int findOrMakeCluster(std::vector<ClusterNode>& n, double x, double y, double tol) {
//    double t2 = tol * tol; for (int i = 0; i < (int)n.size(); ++i)if (dist2(n[i].x, n[i].y, x, y) <= t2)return i; n.push_back({ x,y,{} }); return n.size() - 1;
//}
//static void buildAdjacencyBySharedEndpoints(const std::vector<Segment>& s, double tol, std::vector<std::vector<int>>& adj) {
//    adj.assign(s.size(), {}); std::unordered_map<std::string, std::vector<int>>L; for (int i = 0; i < (int)s.size(); ++i)L[s[i].layer].push_back(i);
//    for (auto& pair : L) {
//        auto& idx = pair.second; std::vector<ClusterNode>nodes;
//        for (int si : idx) {
//            int c1 = findOrMakeCluster(nodes, s[si].startX, s[si].startY, tol); int c2 = findOrMakeCluster(nodes, s[si].endX, s[si].endY, tol);
//            nodes[c1].segs.push_back(si); nodes[c2].segs.push_back(si);
//        }
//        for (auto& nd : nodes)for (size_t a = 0; a < nd.segs.size(); ++a)for (size_t b = a + 1; b < nd.segs.size(); ++b) {
//            int u = nd.segs[a], v = nd.segs[b]; if (u == v)continue; adj[u].push_back(v); adj[v].push_back(u);
//        }
//    }
//}
//static std::vector<std::vector<int>> connectedComponents(const std::vector<std::vector<int>>& adj) {
//    int n = adj.size(); std::vector<int>vis(n, 0); std::vector<std::vector<int>>comps;
//    for (int i = 0; i < n; ++i)if (!vis[i]) {
//        std::vector<int>q{ i }; vis[i] = 1; for (size_t k = 0; k < q.size(); ++k)
//            for (int v : adj[q[k]])if (!vis[v]) { vis[v] = 1; q.push_back(v); }comps.push_back(q);
//    }return comps;
//}
//
//// ====================== 替换与输出 ======================
//static std::string replaceNum(const std::string& blk, const char* key, const std::string& rep) {
//    size_t pos = blk.find("("); while (pos != std::string::npos) {
//        size_t j = pos + 1; while (j < blk.size() && isspace((unsigned char)blk[j]))++j;
//        size_t klen = strlen(key); if (j + klen <= blk.size() && blk.compare(j, klen, key) == 0)break; pos = blk.find("(", pos + 1);
//    }
//    if (pos == std::string::npos)return blk; size_t p = pos + 1; while (isspace((unsigned char)blk[p]))++p; p += strlen(key);
//    while (isspace((unsigned char)blk[p]))++p; size_t nb = p; while (nb < blk.size() && (isdigit((unsigned char)blk[nb]) || blk[nb] == '+' || blk[nb] == '-' || blk[nb] == '.' || blk[nb] == 'e' || blk[nb] == 'E'))++nb;
//    return blk.substr(0, p) + rep + blk.substr(nb);
//}
//static std::string rewriteSmall(const std::string& raw, const std::vector<SegmentSpan>& sp, const std::vector<int>& newNet, const std::vector<double>& newW) {
//    std::string out; size_t cur = 0; char buf[64];
//    for (size_t i = 0; i < sp.size(); ++i) {
//        if (cur < sp[i].begin)out.append(raw, cur, sp[i].begin - cur);
//        std::string blk = raw.substr(sp[i].begin, sp[i].end - sp[i].begin);
//        sprintf_s(buf, sizeof(buf), "%d", newNet[i]);
//        sprintf_s(buf, sizeof(buf), "%.6f", newW[i]);
//        out += blk; cur = sp[i].end;
//    }if (cur < raw.size())out.append(raw, cur, raw.size() - cur); return out;
//}
//static void SaveCSV(const std::string& path, const std::vector<Segment>& s, const std::vector<int>& n, const std::vector<double>& w) {
//    std::ofstream f(path); f << "startX,startY,endX,endY,layer,oldNet,newNet,oldWidth,newWidth\n"; f.setf(std::ios::fixed); f.precision(6);
//    for (size_t i = 0; i < s.size(); ++i)f << s[i].startX << "," << s[i].startY << "," << s[i].endX << "," << s[i].endY << "," << s[i].layer << "," << s[i].net << "," << n[i] << "," << s[i].width << "," << w[i] << "\n";
//}
//
//// ====================== 主程序 ======================
//int main() {
//    try {
//        std::string orig = "C:\\Users\\13437\\Desktop\\testcase.kicad_pcb";
//        std::string small = "C:\\Users\\13437\\Desktop\\output.bestSolutionWithMerging.s_10_i_5_b_0_lc_500_to_15000_vo_40_po_2000_tss_0_vss_0.test1(1).kicad_pcb";
//        double tol = 0.05;
//
//        std::vector<Segment> seg7, seg5;
//        std::vector<SegmentSpan> sp7, sp5;
//        std::string raw;
//        LoadSegmentsWithSpans(orig, seg7, sp7, nullptr);
//        LoadSegmentsWithSpans(small, seg5, sp5, &raw);
//
//        std::cout << "[INFO] 原始(7.0) 段数 = " << seg7.size()
//            << " | rerouting 段数 = " << seg5.size() << "\n";
//
//        MatchParams mp;
//        mp.coordTol = tol;
//        mp.requireLayerEqual = true;
//
//        // 初始化匹配结果
//        std::vector<int> newNet(seg5.size(), -1);
//        std::vector<double> newW(seg5.size(), 0);
//        std::vector<std::vector<std::pair<int, double>>> hits(seg5.size());
//
//        // 逐段匹配
//        for (size_t i = 0; i < seg5.size(); ++i)
//            pickOrigNetAndWidthForSmallSegment(seg5[i], seg7, mp, newNet[i], newW[i], &hits[i]);
//
//        // 聚类路径
//        std::vector<std::vector<int>> adj;
//        buildAdjacencyBySharedEndpoints(seg5, tol, adj);
//        auto comps = connectedComponents(adj);
//
//        // 按路径修正 net 与 width
//        for (auto& comp : comps) {
//            std::map<int, int> voteN;
//            for (int si : comp)
//                for (auto& h : hits[si])
//                    voteN[h.first]++;
//
//            if (voteN.empty()) continue;
//
//            int bestNet = -1, bestVotes = -1;
//            for (auto& kv : voteN)
//                if (kv.second > bestVotes) { bestNet = kv.first; bestVotes = kv.second; }
//
//            std::map<double, int> voteW;
//            for (int si : comp)
//                for (auto& h : hits[si])
//                    if (h.first == bestNet) voteW[h.second]++;
//
//            if (!voteW.empty()) {
//                int wv = -1; double bestW = 0;
//                for (auto& kv : voteW)
//                    if (kv.second > wv) { bestW = kv.first; wv = kv.second; }
//                for (int si : comp) {
//                    newNet[si] = bestNet;
//                    newW[si] = bestW;
//                }
//            }
//        }
//
//        // 生成最终 vector
//        std::vector<Segment> updated;
//        updated.reserve(seg5.size());
//        for (size_t i = 0; i < seg5.size(); ++i) {
//            Segment s = seg5[i];
//            if (newNet[i] != -1) s.net = newNet[i];
//            if (newW[i] != 0) s.width = newW[i];
//            updated.push_back(s);
//        }
//
//        std::cout << "\n 共更新 " << updated.size() << " 条 segment。\n";
//        for (size_t i = 0; i < std::min<size_t>(10, updated.size()); ++i) {
//            const auto& s = updated[i];
//            std::cout << "[" << i + 1 << "] (" << s.startX << "," << s.startY
//                << ") -> (" << s.endX << "," << s.endY << ") "
//                << s.layer << " | width=" << s.width
//                << " | net=" << s.net << "\n";
//        }
//
//        // vector<Segment> updated 就是你要的最终结果
//        // 可在此返回或进一步处理
//    }
//    catch (std::exception& e) {
//        std::cerr << "[ERR] " << e.what() << "\n";
//    }
//}
