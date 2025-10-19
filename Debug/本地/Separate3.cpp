
// 坐标约定：像素 1-based，左上为(1,1)，行向下增大、列向右增大。
// 现实坐标：mm，y 向下增大（与屏幕一致）

#include <iostream>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <string>
#include <iomanip>
#include <algorithm>
#include <array>
#include <climits>
#include "separate3.h"
#include <fstream>
#include <sstream>
#include <map>
#include "First_Part12.h"   // KiCadParser 声明
#include "blocked_area_analyzer.h"   // 内含 grid.h
#include <unordered_set>
#include <unordered_map>
#include <cctype>
#include <filesystem>
#include "cost.h" 
#include <deque>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <climits>  // 为了 INT_MIN
#include <map>      // 你后面用了 std::map





//==========================================SegmentExtractor.cpp==========================================
// ====================== 基础结构 =======================================================================
struct Segment {
    double startX = 0, startY = 0;
    double endX = 0, endY = 0;
    double width = 0;
    std::string layer;
    int net = 0;
};

struct Pad2 {
    double x = 0, y = 0;
    std::string layer;
};



struct SegmentSpan { size_t begin = 0, end = 0; };

static std::string ReadWholeFile(const std::string& path) {
    std::ifstream fin(path, std::ios::binary);
    if (!fin) throw std::runtime_error("无法打开: " + path);
    fin.seekg(0, std::ios::end);
    std::string s(fin.tellg(), '\0');
    fin.seekg(0, std::ios::beg);
    fin.read(&s[0], s.size());
    return s;
}
static void WriteWholeFile(const std::string& path, const std::string& c) {
    std::ofstream f(path, std::ios::binary); if (!f) throw std::runtime_error("写出失败"); f << c;
}
static const char* skip_ws(const char* p) { while (*p && isspace((unsigned char)*p))++p; return p; }
static const char* parse_double(const char* p, double& v) { p = skip_ws(p); char* e = nullptr; v = strtod(p, &e); return e == p ? nullptr : e; }
static const char* parse_int(const char* p, int& v) { p = skip_ws(p); char* e = nullptr; v = strtol(p, &e, 10); return e == p ? nullptr : e; }
static size_t findNextSexpHead(const std::string& t, size_t from, const std::string& kw) {
    for (size_t i = from; i < t.size(); ++i)if (t[i] == '(') {
        size_t j = i + 1; while (j < t.size() && isspace((unsigned char)t[j]))++j;
        if (t.compare(j, kw.size(), kw) == 0)return i;
    }return std::string::npos;
}
static const char* find_kw(const std::string& s, const char* kw) {
    const char* t = s.c_str(); for (size_t i = 0; i < s.size(); ++i) {
        if (t[i] != '(')continue; size_t j = i + 1; while (isspace((unsigned char)t[j]))++j;
        const char* pj = t + j, * pk = kw; while (*pk && *pj == *pk) { ++pj; ++pk; }if (!*pk)return pj;
    }return nullptr;
}
static bool parse_layer_name(const char* p, std::string& o) {
    p = skip_ws(p); if (*p == '"') { ++p; const char* q = p; while (*q && *q != '"')++q; o.assign(p, q); }
    else { const char* q = p; while (*q && !isspace((unsigned char)*q) && *q != ')')++q; o.assign(p, q); }return true;
}
static inline void normalizeLayer(std::string& L) { if (L == "Bottom")L = "B.Cu"; else if (L == "Top")L = "F.Cu"; }
static bool ParseSegmentBlock_NoRegex(const std::string& blk, Segment& s) {
    const char* p;
    if ((p = find_kw(blk, "start"))) { p = parse_double(p, s.startX); p = parse_double(p, s.startY); }
    else return false;
    if ((p = find_kw(blk, "end"))) { p = parse_double(p, s.endX); p = parse_double(p, s.endY); }
    else return false;
    if ((p = find_kw(blk, "width"))) { p = parse_double(p, s.width); }
    else return false;
    if ((p = find_kw(blk, "layer"))) { parse_layer_name(p, s.layer); normalizeLayer(s.layer); }
    else return false;
    if ((p = find_kw(blk, "net"))) { parse_int(p, s.net); }
    else return false; return true;
}
static void LoadSegmentsWithSpans(const std::string& pcb, std::vector<Segment>& outS, std::vector<SegmentSpan>& outSp, std::string* raw = nullptr) {
    std::string t = ReadWholeFile(pcb); if (raw)*raw = t; size_t i = 0; while (true) {
        size_t s = findNextSexpHead(t, i, "segment"); if (s == std::string::npos)break;
        int dep = 0; size_t j = s; bool st = false; for (; j < t.size(); ++j) { if (t[j] == '(') { ++dep; st = true; } else if (t[j] == ')' && st && --dep == 0) { ++j; break; } }
        std::string blk = t.substr(s, j - s); Segment seg; if (ParseSegmentBlock_NoRegex(blk, seg)) { outS.push_back(seg); outSp.push_back({ s,j }); }i = j;
    }
}

static const char* find_kw_local(const std::string& s, const char* kw) {
    size_t pos = s.find(kw);
    if (pos == std::string::npos) return nullptr;
    size_t kwlen = std::char_traits<char>::length(kw);
    return s.c_str() + pos + kwlen; // 返回指到关键字之后的位置
}



// 兼容 KiCad 5.0：读取 (module ...) 的 at 与 layer
static bool ParseModuleAT(const std::string& blk, double& ax, double& ay, double& ang, std::string& layer) {
    const char* p; layer = "F.Cu"; ax = ay = 0; ang = 0;

    // (layer <name>)
    if ((p = find_kw_local(blk, "(layer"))) {
        parse_layer_name(p, layer);
        normalizeLayer(layer);
    }

    // (at x y [angle])
    if (!(p = find_kw_local(blk, "(at"))) return false;
    if (!(p = parse_double(p, ax))) return false;
    if (!(p = parse_double(p, ay))) return false;
    double tmp;
    if ((p = parse_double(p, tmp))) ang = tmp;  // 角度可省略
    return true;
}

// 兼容 KiCad 5.0：读取 (pad ...) 的局部 at（相对 module）
static bool ParsePadLocalAT(const std::string& blk, double& dx, double& dy) {
    const char* p; dx = dy = 0;
    if (!(p = find_kw_local(blk, "(at"))) return false;
    if (!(p = parse_double(p, dx))) return false;
    if (!(p = parse_double(p, dy))) return false;
    return true;
}

// 兼容 KiCad 5.0：把 pad 局部坐标旋转/平移成板上绝对坐标（含层名规范化）
static void LoadPadsAbsolute(const std::string& pcb, std::vector<Pad2>& pads) {
    pads.clear();
    std::string t = ReadWholeFile(pcb);

    size_t i = 0;
    while (true) {
        // 找 (module ...)
        size_t s = findNextSexpHead(t, i, "module");
        if (s == std::string::npos) break;

        // 括号匹配到 module 块末尾
        int dep = 0; size_t j = s; bool st = false;
        for (; j < t.size(); ++j) {
            if (t[j] == '(') { ++dep; st = true; }
            else if (t[j] == ')' && st && --dep == 0) { ++j; break; }
        }
        std::string mblk = t.substr(s, j - s);

        double ax, ay, ang; std::string mLayer;
        if (!ParseModuleAT(mblk, ax, ay, ang, mLayer)) { i = j; continue; }
        double th = ang * M_PI / 180.0;

        // 在 module 内找 (pad ...)
        size_t k = 0;
        while (true) {
            size_t ps = std::string::npos;
            for (size_t u = k; u < mblk.size(); ++u) {
                if (mblk[u] == '(') {
                    size_t v = u + 1;
                    while (v < mblk.size() && isspace((unsigned char)mblk[v])) ++v;
                    if (mblk.compare(v, 3, "pad") == 0) { ps = u; break; }
                }
            }
            if (ps == std::string::npos) break;

            int d2 = 0; size_t pj = ps; bool st2 = false;
            for (; pj < mblk.size(); ++pj) {
                if (mblk[pj] == '(') { ++d2; st2 = true; }
                else if (mblk[pj] == ')' && st2 && --d2 == 0) { ++pj; break; }
            }
            std::string pblk = mblk.substr(ps, pj - ps);

            double dx = 0, dy = 0;
            if (ParsePadLocalAT(pblk, dx, dy)) {
                // 局部 → 绝对
                double xg = ax + std::cos(th) * dx - std::sin(th) * dy;
                double yg = ay + std::sin(th) * dx + std::cos(th) * dy;
                pads.push_back(Pad2{ xg, yg, mLayer });
            }
            k = pj;
        }
        i = j;
    }
}


// ====================== 几何计算 ======================
static inline double sq(double x) { return x * x; }
static inline double dist2(double x1, double y1, double x2, double y2) { return sq(x1 - x2) + sq(y1 - y2); }
static double pointSegmentDist2(double px, double py, double ax, double ay, double bx, double by, bool& on) {
    double vx = bx - ax, vy = by - ay, wx = px - ax, wy = py - ay, vv = vx * vx + vy * vy, t = (vv > 0) ? (wx * vx + wy * vy) / vv : 0; on = true;
    if (t < 0) { on = false; return dist2(px, py, ax, ay); }if (t > 1) { on = false; return dist2(px, py, bx, by); }double cx = ax + t * vx, cy = ay + t * vy; return dist2(px, py, cx, cy);
}
struct MatchParams { double coordTol = 0.05; bool requireLayerEqual = true; };
static bool endpointTouchesOrig(const Segment& s, double px, double py, const Segment& t, const MatchParams& mp) {
    if (mp.requireLayerEqual && s.layer != t.layer)return false; bool on; double d2 = pointSegmentDist2(px, py, t.startX, t.startY, t.endX, t.endY, on);
    return on && d2 <= mp.coordTol * mp.coordTol;
}



//==========================一个“用 pad 做种子”的小函数===============================
// 以 pad 为种子给 small 段预赋 layer/net/width
static void SeedSmallFromPads(
    const std::vector<Pad2>& pads,
    const std::vector<Segment>& orig,   // 原板段(seg7)
    const std::vector<Segment>& small,  // 小板段(seg5)
    double PAD_TOL, double ATTACH_TOL,
    std::vector<int>& newNet,
    std::vector<double>& newW,
    std::vector<std::string>& newLayer)
{
    // 1) 先为每个 pad 找到所落的“原板段”的 (net,width,layer)
    struct PadAttr { int net = -1; double w = 0; std::string layer; bool ok = false; };
    std::vector<PadAttr> pinfo(pads.size());

    for (size_t i = 0; i < pads.size(); ++i) {
        const auto& p = pads[i];
        double best = 1e100; int pick = -1;
        for (int j = 0; j < (int)orig.size(); ++j) {
            if (orig[j].layer != p.layer) continue;
            bool on = false;
            double d2 = pointSegmentDist2(p.x, p.y, orig[j].startX, orig[j].startY, orig[j].endX, orig[j].endY, on);
            if (on && d2 <= PAD_TOL * PAD_TOL && d2 < best) {
                best = d2; pick = j;
            }
        }
        if (pick >= 0) {
            pinfo[i].net = orig[pick].net;
            pinfo[i].w = orig[pick].width;
            pinfo[i].layer = orig[pick].layer;
            pinfo[i].ok = true;
        }
    }

    // 2) 用 pad 给“贴靠到该 pad 的小段端点”做预赋值（同层 + 端点距 pad <= ATTACH_TOL）
    for (size_t i = 0; i < pads.size(); ++i) {
        if (!pinfo[i].ok) continue;
        for (int s = 0; s < (int)small.size(); ++s) {
            if (small[s].layer != pinfo[i].layer) continue;
            if (newNet[s] != INT_MIN && newW[s] >= 0) continue; // 已有更早的赋值就不覆盖
            // 任一端点贴靠 pad 坐标
            if (std::hypot(small[s].startX - pads[i].x, small[s].startY - pads[i].y) <= ATTACH_TOL ||
                std::hypot(small[s].endX - pads[i].x, small[s].endY - pads[i].y) <= ATTACH_TOL)
            {
                newNet[s] = pinfo[i].net;
                newW[s] = pinfo[i].w;
                newLayer[s] = pinfo[i].layer;
            }
        }
    }
}


// ====================== 匹配与聚类 ======================
static void collectHits(const Segment& s, const std::vector<Segment>& orig, const MatchParams& mp, std::vector<std::pair<int, double>>& out) {
    out.clear(); for (auto& t : orig) {
        if (mp.requireLayerEqual && s.layer != t.layer)continue;
        if (endpointTouchesOrig(s, s.startX, s.startY, t, mp))out.emplace_back(t.net, t.width);
        if (endpointTouchesOrig(s, s.endX, s.endY, t, mp))out.emplace_back(t.net, t.width);
    }
}






static void pickOrigNetAndWidthForSmallSegment(const Segment& s, const std::vector<Segment>& orig, const MatchParams& mp,
    int& net, double& w, std::vector<std::pair<int, double>>* outHits = nullptr) {
    std::vector<std::pair<int, double>> hits; collectHits(s, orig, mp, hits); if (outHits)*outHits = hits;
    if (hits.empty()) { net = -1; w = s.width; return; }
    std::map<int, int>vNet; for (auto& h : hits)vNet[h.first]++; net = -1; int bv = -1; for (auto& kv : vNet)
        if (kv.second > bv || (kv.second == bv && kv.first < net)) { bv = kv.second; net = kv.first; }
    std::map<double, int>vW; for (auto& h : hits)if (h.first == net)vW[h.second]++; int wv = -1; w = 0; for (auto& kv : vW)
        if (kv.second > wv || (kv.second == wv && kv.first < w)) { wv = kv.second; w = kv.first; }
}

struct ClusterNode { double x = 0, y = 0; std::vector<int>segs; };
static int findOrMakeCluster(std::vector<ClusterNode>& n, double x, double y, double tol) {
    double t2 = tol * tol; for (int i = 0; i < (int)n.size(); ++i)if (dist2(n[i].x, n[i].y, x, y) <= t2)return i; n.push_back({ x,y,{} }); return n.size() - 1;
}
static void buildAdjacencyBySharedEndpoints(const std::vector<Segment>& s, double tol, std::vector<std::vector<int>>& adj) {
    adj.assign(s.size(), {}); std::unordered_map<std::string, std::vector<int>>L; for (int i = 0; i < (int)s.size(); ++i)L[s[i].layer].push_back(i);
    for (auto& pair : L) {
        auto& idx = pair.second; std::vector<ClusterNode>nodes;
        for (int si : idx) {
            int c1 = findOrMakeCluster(nodes, s[si].startX, s[si].startY, tol); int c2 = findOrMakeCluster(nodes, s[si].endX, s[si].endY, tol);
            nodes[c1].segs.push_back(si); nodes[c2].segs.push_back(si);
        }
        for (auto& nd : nodes)for (size_t a = 0; a < nd.segs.size(); ++a)for (size_t b = a + 1; b < nd.segs.size(); ++b) {
            int u = nd.segs[a], v = nd.segs[b]; if (u == v)continue; adj[u].push_back(v); adj[v].push_back(u);
        }
    }
}
static std::vector<std::vector<int>> connectedComponents(const std::vector<std::vector<int>>& adj) {
    int n = adj.size(); std::vector<int>vis(n, 0); std::vector<std::vector<int>>comps;
    for (int i = 0; i < n; ++i)if (!vis[i]) {
        std::vector<int>q{ i }; vis[i] = 1; for (size_t k = 0; k < q.size(); ++k)
            for (int v : adj[q[k]])if (!vis[v]) { vis[v] = 1; q.push_back(v); }comps.push_back(q);
    }return comps;
}

// ====================== 替换与输出 ======================
// ====================== 替换与输出 ======================
static std::string replaceNum(const std::string& blk, const char* key, const std::string& rep) {
    size_t pos = blk.find("(");
    while (pos != std::string::npos) {
        size_t j = pos + 1;
        while (j < blk.size() && isspace((unsigned char)blk[j])) ++j;
        size_t klen = strlen(key);
        if (j + klen <= blk.size() && blk.compare(j, klen, key) == 0) break;
        pos = blk.find("(", pos + 1);
    }
    if (pos == std::string::npos) return blk;
    size_t p = pos + 1;
    while (isspace((unsigned char)blk[p])) ++p;
    p += strlen(key);
    while (isspace((unsigned char)blk[p])) ++p;

    size_t nb = p;
    while (nb < blk.size() &&
        (isdigit((unsigned char)blk[nb]) || blk[nb] == '+' || blk[nb] == '-' ||
            blk[nb] == '.' || blk[nb] == 'e' || blk[nb] == 'E')) ++nb;

    return blk.substr(0, p) + rep + blk.substr(nb);
}

static std::string rewriteSmall(const std::string& raw,
    const std::vector<SegmentSpan>& sp,
    const std::vector<int>& newNet,
    const std::vector<double>& newW)
{
    std::string out;
    size_t cur = 0;
    char buf[64];

    for (size_t i = 0; i < sp.size(); ++i) {
        if (cur < sp[i].begin) out.append(raw, cur, sp[i].begin - cur);

        std::string blk = raw.substr(sp[i].begin, sp[i].end - sp[i].begin);

        // 仅当提供了新值时才替换；否则保留原文本
        if (i < newNet.size() && newNet[i] != -1) {
            std::snprintf(buf, sizeof(buf), "%d", newNet[i]);
            blk = replaceNum(blk, "net", buf);
        }
        if (i < newW.size() && newW[i] > 0) {
            std::snprintf(buf, sizeof(buf), "%.6f", newW[i]);
            blk = replaceNum(blk, "width", buf);
        }

        out += blk;
        cur = sp[i].end;
    }
    if (cur < raw.size()) out.append(raw, cur, raw.size() - cur);
    return out;
}

static void SaveCSV(const std::string& path, const std::vector<Segment>& s, const std::vector<int>& n, const std::vector<double>& w) {
    std::ofstream f(path); f << "startX,startY,endX,endY,layer,oldNet,newNet,oldWidth,newWidth\n"; f.setf(std::ios::fixed); f.precision(6);
    for (size_t i = 0; i < s.size(); ++i)f << s[i].startX << "," << s[i].startY << "," << s[i].endX << "," << s[i].endY << "," << s[i].layer << "," << s[i].net << "," << n[i] << "," << s[i].width << "," << w[i] << "\n";
}

//========================================================================================================

// 打印从 KiCad 5.0 解析出的、需要加入的新线段（updated）
static void PrintSegmentsToAddFromV5(const std::vector<Segment>& updated) {
    std::cout << "[V5 需新增线段] 共 " << updated.size() << " 条\n";
    for (size_t i = 0; i < updated.size(); ++i) {
        const auto& s = updated[i];
        std::cout << "  #" << i
            << " start(" << s.startX << ", " << s.startY << ")"
            << " end(" << s.endX << ", " << s.endY << ")"
            << " width=" << s.width
            << " layer=" << s.layer
            << " net=" << s.net
            << "\n";
    }
}




// —— 小工具：从 [[SEID:x]] 里解析出 x；失败返回 -1
static int ExtractSEID(const std::shared_ptr<Node>& segment) {
    if (!segment || segment->name != "segment" || segment->parameters.empty())
        return -1;
    const std::string& p0 = segment->parameters.front(); // 期望 "[[SEID:x]]"
    const std::string prefix = "[[SEID:";
    const std::string suffix = "]]";
    if (p0.rfind(prefix, 0) != 0) return -1;
    size_t colon = p0.find(':');
    size_t end = p0.rfind(suffix);
    if (colon == std::string::npos || end == std::string::npos || end <= colon + 1) return -1;
    try {
        return std::stoi(p0.substr(colon + 1, end - colon - 1));
    }
    catch (...) { return -1; }
}


//--------------------------------------------删除end或者start在矩形1里面的所有线段-------------------------------------
// —— 收集 SEID：只看“第一个矩形”（像素域扩张 → 现实坐标外沿），含边界
std::vector<int>
CollectSegmentsTouchingFirstRect(const Grid& grid, int layerId, int blockedAreaId, int expandPixels) {
    std::vector<int> out;

    // 第一个矩形：你第一步“像素拓展”后得到的现实外沿矩形
    // 现有函数 GetExpandedPCBRealCorners(layerId, blockedAreaId, expandPixels)
    // 返回顺序 LT, RT, RB, LB，我们转成矩形盒子
    auto corners = GetExpandedPCBRealCorners(grid, layerId, blockedAreaId, expandPixels);
    if (corners.size() != 4) return out;
    RealRectMM firstRect = {
        /*x_left*/  corners[0].first,
        /*y_top*/   corners[0].second,
        /*x_right*/ corners[1].first,
        /*y_bot*/   corners[2].second
    };

    const std::string targetLayer = getLayerNameById(grid, layerId);

    std::unordered_set<int> uniq;//创建一个无序的set用来存储segment的序号：set的名字是uniq
    uniq.reserve(grid.segments.size());

    for (const auto& seg : grid.segments) {
        double x1, y1, x2, y2;
        std::string layer, net;
        if (!GetSegmentInfo(seg, x1, y1, x2, y2, layer, net)) continue;
        if (layer != targetLayer) continue;

        // “在里面或边上” → IsPointInsideRect() 本身就是闭区间判断（>= && <=）
        if (IsPointInsideRect(x1, y1, firstRect) || IsPointInsideRect(x2, y2, firstRect)) {
            int id = ExtractSEID(seg);
            if (id > 0) uniq.insert(id);//将ID存入unordered_set（真正的存储）
        }
    }

    out.assign(uniq.begin(), uniq.end());// 从set复制到vector（数据转移）
    std::sort(out.begin(), out.end());
    return out;
}

// —— 批量删除：用 First_Part12.cpp 里的删除接口
size_t
DeleteSegmentsTouchingFirstRect(KiCadParser& parser, const Grid& grid, int layerId, int blockedAreaId, int expandPixels) {
    auto ids = CollectSegmentsTouchingFirstRect(grid, layerId, blockedAreaId, expandPixels);
    size_t removedTotal = 0;
    for (int id : ids) {
        // KiCadParser::removeSegmentById 会在内部重建缓存（见 First_Part12.cpp）
        removedTotal += parser.removeSegmentById(id);
    }
    return removedTotal;
}

//--------------------------------------------删除end或者start在矩形1里面的所有线段-------------------------------------








// ========== 像素中心 <-> 现实坐标 ==========
std::pair<double, double>
PixelCenterPx1ToRealMM(const Grid& grid, int row1, int col1)
{
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    row1 = std::max(1, std::min(row1, grid.height));
    col1 = std::max(1, std::min(col1, grid.width));
    const double s = (double)grid.inputScale;
    const double x = grid.min_x + ((double)col1 - 0.5) / s;
    const double y = grid.min_y + ((double)row1 - 0.5) / s; // 向下增大
    return { x, y };
}

std::pair<int, int>
RealMMToPixelPx1(const Grid& grid, double x_mm, double y_mm)
{
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    const double s = (double)grid.inputScale;
    int col1 = (int)std::floor((x_mm - grid.min_x) * s) + 1;
    int row1 = (int)std::floor((y_mm - grid.min_y) * s) + 1;
    col1 = std::max(1, std::min(col1, grid.width));
    row1 = std::max(1, std::min(row1, grid.height));
    return { row1, col1 };
}
// ========== 像素边界 -> 现实外沿矩形 ==========
//给出边界像素坐标，输出外延现实坐标
RealRectMM PixelBoxEdgesPx1ToRealRect(const Grid& grid,
    int row_min, int col_min,
    int row_max, int col_max)
{
    if (grid.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    const double s = (double)grid.inputScale;
    RealRectMM rr;
    rr.x_left = grid.min_x + ((double)col_min - 1.0) / s;
    rr.x_right = grid.min_x + ((double)col_max) / s;
    rr.y_top = grid.min_y + ((double)row_min - 1.0) / s;
    rr.y_bot = grid.min_y + ((double)row_max) / s;
    return rr;
}

// ========== 获取阻塞像素（1-based） ==========
std::vector<std::pair<int, int>>
GetBlockedAreaPixels_1Based(const Grid& grid, int layer, int blockedAreaID)
{
    if (layer < 0 || layer >= grid.Layers)
        throw std::runtime_error("层号超出范围");

    LayerBlockedAnalysis la = BlockedAreaAnalyzer::analyzeLayerBlockedAreas(grid, layer);

    const BlockedAreaInfo* target = nullptr;
    for (const auto& area : la.blockedAreas) {
        if (area.areaID == blockedAreaID) { target = &area; break; }
    }
    if (!target) return {};

    std::vector<std::pair<int, int>> out;
    out.reserve(target->cells.size());
    for (const auto& rc0 : target->cells) {
        out.emplace_back(rc0.first + 1, rc0.second + 1); // 0-based -> 1-based
    }
    return out;
}


// ========== 计算像素边界框（1-based） ==========
BBoxPx1 ComputeBBoxFromPixels1(const std::vector<std::pair<int, int>>& pxs1)
{
    if (pxs1.empty()) throw std::runtime_error("ComputeBBoxFromPixels1: 输入像素为空");
    int rmin = INT_MAX, cmin = INT_MAX, rmax = INT_MIN, cmax = INT_MIN;
    for (auto [r1, c1] : pxs1) {
        rmin = std::min(rmin, r1);
        cmin = std::min(cmin, c1);
        rmax = std::max(rmax, r1);
        cmax = std::max(cmax, c1);
    }
    return { rmin, cmin, rmax, cmax };
}

// 像素四角（LB, RB, RT, LT）
std::array<std::pair<int, int>, 4> BBoxCorners_LB_RB_RT_LT(const BBoxPx1& b)
{
    return {
        std::pair<int,int>{ b.row_max, b.col_min }, // LB
        std::pair<int,int>{ b.row_max, b.col_max }, // RB
        std::pair<int,int>{ b.row_min, b.col_max }, // RT
        std::pair<int,int>{ b.row_min, b.col_min }  // LT
    };
}





// 现实四角（严格 LT, LB, RT, RB，用于切板/绘制边框）
std::array<std::pair<double, double>, 4>
RealCorners_LT_LB_RT_RB(const Grid& grid, const BBoxPx1& b)
{
    const RealRectMM rr = PixelBoxEdgesPx1ToRealRect(grid, b.row_min, b.col_min, b.row_max, b.col_max);



    return {
        std::pair<double,double>{ rr.x_left , rr.y_top }, // LT
        std::pair<double,double>{ rr.x_left , rr.y_bot }, // LB
        std::pair<double,double>{ rr.x_right, rr.y_top }, // RT
        std::pair<double,double>{ rr.x_right, rr.y_bot }  // RB
    };
}

//// 现实四角（使用像素中心坐标，严格 LT, LB, RT, RB，用于切板/绘制边框）
//std::array<std::pair<double, double>, 4>
//RealCorners_LT_LB_RT_RB(const Grid& grid, const BBoxPx1& b)
//{
//    // 直接使用第一个函数计算四个角的像素中心坐标
//    return {
//        PixelCenterPx1ToRealMM(grid, b.row_min, b.col_min), // 左上角 (LT)
//        PixelCenterPx1ToRealMM(grid, b.row_max, b.col_min), // 左下角 (LB)
//        PixelCenterPx1ToRealMM(grid, b.row_min, b.col_max), // 右上角 (RT)
//        PixelCenterPx1ToRealMM(grid, b.row_max, b.col_max)  // 右下角 (RB)
//    };
//}



// 逻辑：希望左右各扩 m、上下各扩 m；若某侧不够扩，把缺的补到对侧。
BBoxPx1 ExpandBBoxByMPixelsCompensated(const BBoxPx1& b, int m,
    const Grid& g, ClipReport* rep)
{
    const int W = g.width, H = g.height;

    // 该区域到边界的“可扩”空间
    const int avail_left = b.col_min - 1;
    const int avail_right = W - b.col_max;
    const int avail_top = b.row_min - 1;
    const int avail_bottom = H - b.row_max;

    const int left_def = std::max(0, m - avail_left);
    const int right_def = std::max(0, m - avail_right);
    const int top_def = std::max(0, m - avail_top);
    const int bot_def = std::max(0, m - avail_bottom);

    if (rep) {
        rep->left = (left_def > 0);
        rep->right = (right_def > 0);
        rep->top = (top_def > 0);
        rep->bottom = (bot_def > 0);
        rep->left_deficit_px = left_def;
        rep->right_deficit_px = right_def;
        rep->top_deficit_px = top_def;
        rep->bottom_deficit_px = bot_def;
    }

    // 对侧补偿：希望总扩张量 ≈ 2m
    const int extend_left = m + right_def; // 右侧不够，把缺的补到左侧
    const int extend_right = m + left_def;  // 左侧不够，把缺的补到右侧
    const int extend_top = m + bot_def;   // 下侧不够，补到上侧
    const int extend_bottom = m + top_def;   // 上侧不够，补到下侧

    BBoxPx1 e;
    e.col_min = std::max(1, b.col_min - extend_left);
    e.col_max = std::min(W, b.col_max + extend_right);
    e.row_min = std::max(1, b.row_min - extend_top);
    e.row_max = std::min(H, b.row_max + extend_bottom);

    // 防翻转（极端窄板+超大 m）
    if (e.col_min > e.col_max) { e.col_min = 1; e.col_max = W; }
    if (e.row_min > e.row_max) { e.row_min = 1; e.row_max = H; }

    return e;
}


// 根据层编号获取层名称
std::string getLayerNameById(const Grid& grid, int layerId) {
    if (layerId >= 0 && layerId < grid.layerName.size()) {
        return grid.layerName[layerId]; // 返回对应层编号的层名称
    }
    else {
        return "无效的层编号"; // 如果层编号超出范围，返回错误信息
    }
}


// 根据层名称获取层编号
int getLayerIdByName(const Grid& grid, const std::string& layerName) {
    for (int i = 0; i < grid.layerName.size(); i++) {
        if (grid.layerName[i] == layerName) {
            return i; // 返回对应的层编号
        }
    }
    return -1; // 如果找不到，返回-1
}






static inline RealRectMM ExpandRect_HalfPlusK(const RealRectMM& in, const Grid& grid, int kWholePixels) {
    if (grid.inputScale <= 0) throw std::runtime_error("inputScale 必须为正");
    const double px_mm = 1.0 / static_cast<double>(grid.inputScale); // 1px = 1/inputScale mm
    const double d = (static_cast<double>(kWholePixels)) * px_mm;
    RealRectMM r = in;
    r.x_left -= d;  r.x_right += d;
    r.y_top -= d;  r.y_bot += d;
    return r;
}

static inline RealRectMM ClipRectToBoard(const RealRectMM& in,
    double boardMinX, double boardMinY,
    double boardMaxX, double boardMaxY) {
    RealRectMM r = in;
    if (r.x_left < boardMinX) r.x_left = boardMinX;
    if (r.x_right > boardMaxX) r.x_right = boardMaxX;
    if (r.y_top < boardMinY) r.y_top = boardMinY;
    if (r.y_bot > boardMaxY) r.y_bot = boardMaxY;
    // 极端防翻转
    if (r.x_left > r.x_right) { r.x_left = boardMinX; r.x_right = boardMaxX; }
    if (r.y_top > r.y_bot) { r.y_top = boardMinY; r.y_bot = boardMaxY; }
    return r;
}

// 方案A主函数（建议用这个）：像素域扩张 -> 现实域 0.5px+kpx -> 裁剪板框 -> 返回 LT,RT,RB,LB
std::vector<std::pair<double, double>>
GetExpandedPCBCornersSnapAndClip(const Grid& grid, int layer, int blockedAreaID,
    int expandPixels, int kWholePixels,
    double boardMinX, double boardMinY,
    double boardMaxX, double boardMaxY)
{
    // 1) 先得到“像素域扩张后的现实四角”（LT,RT,RB,LB）
    auto baseCorners = GetExpandedPCBRealCorners(grid, layer, blockedAreaID, expandPixels);
    if (baseCorners.size() != 4) {
        throw std::runtime_error("GetExpandedPCBCornersSnapAndClip: 角点数 != 4");
    }

    // 2) 现实域对称外扩 0.5px + k×1px
    RealRectMM baseRect = CornersToRect_LT_RT_RB_LB(baseCorners);
    RealRectMM grown = ExpandRect_HalfPlusK(baseRect, grid, kWholePixels);

    // 3) 裁剪到真实板框
    RealRectMM clipped = ClipRectToBoard(grown, boardMinX, boardMinY, boardMaxX, boardMaxY);

    // 4) 回四角（LT,RT,RB,LB）
    return RectToCorners_LT_RT_RB_LB(clipped);
}


// ========== 获取扩展后小PCB的四个现实角点坐标 ==========
std::vector<std::pair<double, double>>
GetExpandedPCBRealCorners(const Grid& grid, int layer, int blockedAreaID, int expandPixels)
{



    // 1. 获取阻塞区域像素
    auto blockedPixels = GetBlockedAreaPixels_1Based(grid, layer, blockedAreaID);
    if (blockedPixels.empty()) {
        throw std::runtime_error("未找到指定的阻塞区域");
    }

    // 2. 计算原始包围盒
    BBoxPx1 originalBBox = ComputeBBoxFromPixels1(blockedPixels);

    // 3. 扩展包围盒
    ClipReport clipReport;
    BBoxPx1 expandedBBox = ExpandBBoxByMPixelsCompensated(originalBBox, expandPixels, grid, &clipReport);

    // 4. 获取现实坐标的四个角点（顺时针顺序：LT, RT, RB, LB）
    auto realCorners = RealCorners_LT_LB_RT_RB(grid, expandedBBox);

    // 5. 重新排列为顺时针顺序：LT -> RT -> RB -> LB
    std::vector<std::pair<double, double>> result;
    result.reserve(4);

    // LT (左上)
    result.push_back(realCorners[0]);
    // RT (右上)  
    result.push_back(realCorners[2]);
    // RB (右下)
    result.push_back(realCorners[3]);
    // LB (左下)
    result.push_back(realCorners[1]);

    return result;
}

// 使用示例：
void exampleUsage(const Grid& grid) {
    try {
        int layer = 0;          // 层号
        int blockedAreaID = 1;  // 阻塞区域ID
        int expandPixels = 5;   // 扩展像素数

        auto corners = GetExpandedPCBRealCorners(grid, layer, blockedAreaID, expandPixels);

        // 打印四个角点坐标
        std::cout << "小PCB的四个角点坐标（顺时针顺序）：" << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        for (size_t i = 0; i < corners.size(); ++i) {
            std::cout << "角点 " << i + 1 << ": ("
                << corners[i].first << " mm, "
                << corners[i].second << " mm)" << std::endl;
        }

        // 如果需要直接使用这个vector
        // corners 现在包含了四个点的坐标，按顺时针顺序：LT->RT->RB->LB

    }
    catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
    }
}



// ========== 计算线段与矩形边界的交点 ==========
std::vector<std::pair<double, double>>
ComputeLineRectIntersections(double x1, double y1, double x2, double y2, const RealRectMM& rect)
{
    std::vector<std::pair<double, double>> intersections;

    // 检查与四条边的交点
    // 左边界 (x = rect.x_left)
    if (x1 != x2) {
        double t = (rect.x_left - x1) / (x2 - x1);
        if (t >= 0 && t <= 1) {
            double y = y1 + t * (y2 - y1);
            if (y >= rect.y_top && y <= rect.y_bot) {
                intersections.emplace_back(rect.x_left, y);
            }
        }
    }

    // 右边界 (x = rect.x_right)
    if (x1 != x2) {
        double t = (rect.x_right - x1) / (x2 - x1);
        if (t >= 0 && t <= 1) {
            double y = y1 + t * (y2 - y1);
            if (y >= rect.y_top && y <= rect.y_bot) {
                intersections.emplace_back(rect.x_right, y);
            }
        }
    }

    // 上边界 (y = rect.y_top)
    if (y1 != y2) {
        double t = (rect.y_top - y1) / (y2 - y1);
        if (t >= 0 && t <= 1) {
            double x = x1 + t * (x2 - x1);
            if (x >= rect.x_left && x <= rect.x_right) {
                intersections.emplace_back(x, rect.y_top);
            }
        }
    }

    // 下边界 (y = rect.y_bot)
    if (y1 != y2) {
        double t = (rect.y_bot - y1) / (y2 - y1);
        if (t >= 0 && t <= 1) {
            double x = x1 + t * (x2 - x1);
            if (x >= rect.x_left && x <= rect.x_right) {
                intersections.emplace_back(x, rect.y_bot);
            }
        }
    }

    return intersections;
}

// ========== 判断点是否在矩形内部 ==========
bool IsPointInsideRect(double x, double y, const RealRectMM& rect)
{
    return (x >= rect.x_left && x <= rect.x_right &&
        y >= rect.y_top && y <= rect.y_bot);
}

// ========== 修改后的获取线段信息函数 ==========
bool GetSegmentInfo(const std::shared_ptr<Node>& segment,
    double& x1, double& y1, double& x2, double& y2,
    std::string& layer, std::string& net)
{
    x1 = y1 = x2 = y2 = 0.0;
    layer.clear();
    net.clear();

    bool hasStart = false, hasEnd = false, hasLayer = false, hasNet = false;

    for (const auto& child : segment->children) {
        if (child->name == "start") {
            if (child->parameters.size() >= 2) {
                x1 = std::stod(child->parameters[0]);
                y1 = std::stod(child->parameters[1]);
                hasStart = true;
            }
        }
        else if (child->name == "end") {
            if (child->parameters.size() >= 2) {
                x2 = std::stod(child->parameters[0]);
                y2 = std::stod(child->parameters[1]);
                hasEnd = true;
            }
        }
        else if (child->name == "layer") {
            if (!child->parameters.empty()) {
                layer = child->parameters[0];
                hasLayer = true;
            }
        }
        else if (child->name == "net") {
            if (!child->parameters.empty()) {
                net = child->parameters[0];
                hasNet = true;
            }
        }
    }

    return (hasStart && hasEnd && hasLayer && hasNet);
}

// ========== 修改后的判断连接点函数 ==========
bool IsConnectionPoint(double x, double y, const std::vector<std::shared_ptr<Node>>& segments,
    const std::string& targetLayer, double tolerance)
{
    int connectionCount = 0;

    for (const auto& segment : segments) {
        double x1, y1, x2, y2;
        std::string layer, net;

        if (!GetSegmentInfo(segment, x1, y1, x2, y2, layer, net)) {
            continue;
        }

        // 检查是否在目标层
        if (layer != targetLayer) {
            continue;
        }

        // 检查点是否与线段起点匹配
        if (std::abs(x1 - x) < tolerance && std::abs(y1 - y) < tolerance) {
            connectionCount++;
        }

        // 检查点是否与线段终点匹配
        if (std::abs(x2 - x) < tolerance && std::abs(y2 - y) < tolerance) {
            connectionCount++;
        }

        // 如果连接数超过1，说明是连接点
        if (connectionCount > 1) {
            return true;
        }
    }

    return false;
}





std::vector<PointWithNet>
FindIntersectionsAndEndpointsInSmallPCB(const Grid& grid, int layerId, int blockedAreaId, int expandPixels)
{
    std::vector<PointWithNet> result;

    try {
        // 1. 获取小PCB的边界框
        auto corners = GetExpandedPCBRealCorners(grid, layerId, blockedAreaId, expandPixels);
        RealRectMM pcbBox;
        pcbBox.x_left = corners[0].first;
        pcbBox.y_top = corners[0].second;
        pcbBox.x_right = corners[1].first;
        pcbBox.y_bot = corners[2].second;

        std::cout << "小PCB边界框: 左=" << pcbBox.x_left << ", 上=" << pcbBox.y_top
            << ", 右=" << pcbBox.x_right << ", 下=" << pcbBox.y_bot << std::endl;

        // 2. 获取当前层的名称
        std::string targetLayerName = getLayerNameById(grid, layerId);
        std::cout << "目标层: " << targetLayerName << std::endl;

        // 3. 直接使用 grid.segments
        const auto& segments = grid.segments;
        std::cout << "找到 " << segments.size() << " 个线段" << std::endl;

        // 4. 遍历所有线段
        int processedCount = 0;
        for (const auto& segment : segments) {
            double x1, y1, x2, y2;
            std::string layer, net;

            if (!GetSegmentInfo(segment, x1, y1, x2, y2, layer, net)) {
                continue;
            }

            // 检查是否在目标层
            if (layer != targetLayerName) {
                continue;
            }

            processedCount++;
            std::cout << "处理线段 " << processedCount << ": (" << x1 << "," << y1 << ") -> (" << x2 << "," << y2
                << "), 层: " << layer << ", 网络: " << net << std::endl;

            // 5. 计算线段与小PCB边界的交点
            auto intersections = ComputeLineRectIntersections(x1, y1, x2, y2, pcbBox);

            // 6. 检查线段端点是否在小PCB内部
            bool startInside = IsPointInsideRect(x1, y1, pcbBox);
            bool endInside = IsPointInsideRect(x2, y2, pcbBox);

            std::cout << "  起点在内部: " << (startInside ? "是" : "否")
                << ", 终点在内部: " << (endInside ? "是" : "否")
                << ", 交点数量: " << intersections.size() << std::endl;

            // 7. 收集交点（不检查是否为连接点，因为交点在边界上）
            for (const auto& intersection : intersections) {
                // 交点通常在边界上，不需要检查是否为连接点
                result.emplace_back(intersection.first, intersection.second, net);
                std::cout << "  交点: (" << intersection.first << ", " << intersection.second
                    << "), 网络: " << net << std::endl;
            }

            // 8. 检查终止点，排除连接点
            if (startInside && !endInside) {
                // 起点在内部，终点在外部 → 起点是终止点
                // 检查起点是否为连接点
                if (!IsConnectionPoint(x1, y1, segments, targetLayerName)) {
                    result.emplace_back(x1, y1, net);
                    std::cout << "  终止点(起点在内部): (" << x1 << ", " << y1
                        << "), 网络: " << net << std::endl;
                }
                else {
                    std::cout << "  排除连接点(起点): (" << x1 << ", " << y1
                        << "), 网络: " << net << std::endl;
                }
            }
            else if (!startInside && endInside) {
                // 起点在外部，终点在内部 → 终点是终止点
                // 检查终点是否为连接点
                if (!IsConnectionPoint(x2, y2, segments, targetLayerName)) {
                    result.emplace_back(x2, y2, net);
                    std::cout << "  终止点(终点在内部): (" << x2 << ", " << y2
                        << "), 网络: " << net << std::endl;
                }
                else {
                    std::cout << "  排除连接点(终点): (" << x2 << ", " << y2
                        << "), 网络: " << net << std::endl;
                }
            }
            else if (startInside && endInside) {
                // 线段完全在内部，两个端点都在小PCB内
                // 检查两个端点是否为连接点
                if (!IsConnectionPoint(x1, y1, segments, targetLayerName)) {
                    result.emplace_back(x1, y1, net);
                    std::cout << "  终止点(起点在内部): (" << x1 << ", " << y1
                        << "), 网络: " << net << std::endl;
                }
                else {
                    std::cout << "  排除连接点(起点): (" << x1 << ", " << y1
                        << "), 网络: " << net << std::endl;
                }

                if (!IsConnectionPoint(x2, y2, segments, targetLayerName)) {
                    result.emplace_back(x2, y2, net);
                    std::cout << "  终止点(终点在内部): (" << x2 << ", " << y2
                        << "), 网络: " << net << std::endl;
                }
                else {
                    std::cout << "  排除连接点(终点): (" << x2 << ", " << y2
                        << "), 网络: " << net << std::endl;
                }
            }
        }

        std::cout << "共处理了 " << processedCount << " 个在目标层上的线段" << std::endl;

        // 9. 去重（确保结果中没有重复点）
        auto last = std::unique(result.begin(), result.end(),
            [](const PointWithNet& a, const PointWithNet& b) {
                return std::abs(a.x - b.x) < 0.001 &&
                    std::abs(a.y - b.y) < 0.001 &&
                    a.net == b.net;
            });
        result.erase(last, result.end());

    }
    catch (const std::exception& e) {
        std::cerr << "在FindIntersectionsAndEndpointsInSmallPCB中出错: " << e.what() << std::endl;
    }

    return result;
}


// === 像素大小（mm/px） ===
static inline double PixelSizeMM(const Grid& g) {
    if (g.inputScale <= 0) throw std::runtime_error("grid.inputScale 必须为正");
    return 1.0 / static_cast<double>(g.inputScale);
}

// === 用四角（LT,RT,RB,LB）组装矩形 ===
static inline RealRectMM CornersToRect_LT_RT_RB_LB(const std::vector<std::pair<double, double>>& c) {
    if (c.size() != 4) throw std::runtime_error("CornersToRect_LT_RT_RB_LB: 需要4个角点");
    RealRectMM r;
    // 约定：c[0]=LT, c[1]=RT, c[2]=RB, c[3]=LB
    r.x_left = c[0].first;
    r.y_top = c[0].second;
    r.x_right = c[1].first;
    r.y_bot = c[2].second;
    return r;
}

// === 把矩形转回四角（LT,RT,RB,LB） ===
static inline std::vector<std::pair<double, double>> RectToCorners_LT_RT_RB_LB(const RealRectMM& r) {
    return {
        { r.x_left , r.y_top },   // LT
        { r.x_right, r.y_top },   // RT
        { r.x_right, r.y_bot },   // RB
        { r.x_left , r.y_bot }    // LB
    };
}






// —— 辅助：把“在边界上的现实坐标”稳定地映射到小板内部的像素 ——
// 注意：右边界/下边界要往里收一个极小量，否则会落到外侧像素。
static inline std::pair<int, int> RealToPixel_Inside_LT(const Grid& g, double x, double y) {
    return RealMMToPixelPx1(g, x, y); // 左/上边界直接 floor 落在内侧
}
static inline std::pair<int, int> RealToPixel_Inside_RT(const Grid& g, double x, double y) {
    const double eps = 1.0 / (g.inputScale * 1024.0); // < 1px 的很小量
    return RealMMToPixelPx1(g, x - eps, y);
}
static inline std::pair<int, int> RealToPixel_Inside_RB(const Grid& g, double x, double y) {
    const double eps = 1.0 / (g.inputScale * 1024.0);
    return RealMMToPixelPx1(g, x - eps, y - eps);
}
static inline std::pair<int, int> RealToPixel_Inside_LB(const Grid& g, double x, double y) {
    const double eps = 1.0 / (g.inputScale * 1024.0);
    return RealMMToPixelPx1(g, x, y - eps);
}

// —— 对外1：四个角的“锚点像素”（LT,RT,RB,LB） ——
// 先拿现实角点(LT,RT,RB,LB)，再做“内侧偏置”的像素映射。
std::array<std::pair<int, int>, 4>
GetExpandedPCBCornerAnchorPixels(const Grid& grid,
    int layer, int blockedAreaID,
    int expandPixels)
{
    auto corners = GetExpandedPCBRealCorners(grid, layer, blockedAreaID, expandPixels);
    if (corners.size() != 4) throw std::runtime_error("角点数量应为4");

    // corners[0]=LT, [1]=RT, [2]=RB, [3]=LB （当前实现返回顺序即如此）
    auto lt = RealToPixel_Inside_LT(grid, corners[0].first, corners[0].second);
    auto rt = RealToPixel_Inside_RT(grid, corners[1].first, corners[1].second);
    auto rb = RealToPixel_Inside_RB(grid, corners[2].first, corners[2].second);
    auto lb = RealToPixel_Inside_LB(grid, corners[3].first, corners[3].second);
    return { lt, rt, rb, lb };
}

// —— 辅助：根据“角的方向”生成 N×N 的紫色块，并做 1..H/1..W 裁剪 ——
// dr_sign: 行方向(+1向下, -1向上)；dc_sign: 列方向(+1向右, -1向左)
static inline void AppendCornerPatch(std::vector<std::pair<int, int>>& out,
    int r0, int c0, int N,
    int dr_sign, int dc_sign,
    const Grid& g)
{
    for (int di = 0; di < N; ++di) {
        int r = r0 + dr_sign * di;
        if (r < 1 || r > g.height) continue;
        for (int dj = 0; dj < N; ++dj) {
            int c = c0 + dc_sign * dj;
            if (c < 1 || c > g.width) continue;
            out.emplace_back(r, c);
        }
    }
}

// —— 对外2：四个角的紫色块像素，按 LT→RT→RB→LB 顺序拼接到一个 vector ——
// 例：patchN=3 就是每个角 3×3；若只要一个像素，传 patchN=1 即可。
std::vector<std::pair<int, int>>
GetCornerPixelsPx1(const Grid& grid,
    int layer, int blockedAreaID,
    int expandPixels,
    int patchN)
{
    if (patchN < 1) patchN = 1;

    auto anchors = GetExpandedPCBCornerAnchorPixels(grid, layer, blockedAreaID, expandPixels);
    std::vector<std::pair<int, int>> out;
    out.reserve(4 * patchN * patchN);

    // LT：向“下、右”扩 N×N
    AppendCornerPatch(out, anchors[0].first, anchors[0].second, patchN, +1, +1, grid);
    // RT：向“下、左”扩
    AppendCornerPatch(out, anchors[1].first, anchors[1].second, patchN, +1, -1, grid);
    // RB：向“上、左”扩
    AppendCornerPatch(out, anchors[2].first, anchors[2].second, patchN, -1, -1, grid);
    // LB：向“上、右”扩
    AppendCornerPatch(out, anchors[3].first, anchors[3].second, patchN, -1, +1, grid);

    return out;
}

// ========== 生成完整PCB模板文件 ==========
void GenerateCompletePCBTemplate(
    const std::vector<PointWithNet>& modulePoints,
    const std::vector<std::pair<double, double>>& expandedCorners,
    const std::string& outputFilename)
{
    // 四个module的模板
    std::string template_four_modules = R"((kicad_pcb (version 20171130) (host pcbnew "(5.1.2)-2")

  (general
    (thickness 1.6)
  )

  (page A4)
  (layers
    (0 Top signal)
    (31 Bottom signal)
    (32 B.Adhes user)
    (33 F.Adhes user)
    (34 B.Paste user)
    (35 F.Paste user)
    (36 B.SilkS user)
    (37 F.SilkS user)
    (38 B.Mask user)
    (39 F.Mask user)
    (40 Dwgs.User user)
    (41 Cmts.User user)
    (42 Eco1.User user)
    (43 Eco2.User user)
    (44 Edge.Cuts user)
    (45 Margin user)
    (46 B.CrtYd user)
    (47 F.CrtYd user)
    (48 B.Fab user)
    (49 F.Fab user)
  )

  (setup
    (pad_to_mask_clearance 0.1016)
    (aux_axis_origin -120.62617 153.22756)
    (grid_origin -120.62617 153.22756)
    (pcbplotparams
      (layerselection 0x00010fc_ffffffff)
      (plot_on_all_layers_selection 0x0000000_00000000)
      (disableapertmacros false)
      (usegerberextensions false)
      (usegerberattributes true)
      (usegerberadvancedattributes true)
      (creategerberjobfile true)
      (dashed_line_dash_ratio 12)
      (dashed_line_gap_ratio 3)
      (svgprecision 4)
      (plotframeref false)
      (viasonmask false)
      (mode 1)
      (useauxorigin false)
      (hpglpennumber 1)
      (hpglpenspeed 20)
      (hpglpendiameter 15)
      (dxfpolygonmode true)
      (dxfimperialunits true)
      (dxfusepcbnewfont true)
      (psnegative false)
      (psa4output false)
      (plotreference true)
      (plotvalue true)
      (plotinvisibletext false)
      (sketchpadsonfab false)
      (subtractmaskfromsilk false)
      (outputformat 1)
      (mirror false)
      (drillshape 1)
      (scaleselection 1)
      (outputdirectory "")
    )
  )

  (net 0 "")

  (net 1 Y+)
  (net 2 Y-)
  (net 3 NetR3_2)
  (net 4 W-)
  (net 5 W+)
  (net 6 VIN)
  (net 7 VCC_5V)
  (net 8 TK6)
  (net 9 TK5)
  (net 10 TK4)
  (net 11 TK3)
  (net 12 TK2)
  (net 13 TK1)
  (net 14 SK4)
  (net 15 SK3)
  (net 16 SK2)
  (net 17 SK1)
  (net 18 S_LED7)
  (net 19 S_LED6)
  (net 20 S_LED5)
  (net 21 S_LED4)
  (net 22 S_LED3)
  (net 23 S_LED2)
  (net 24 S_LED1)
  (net 25 PWM2)
  (net 26 PWM1)
  (net 27 NetR19_2)
  (net 28 NetR16_2)
  (net 29 NetR15_2)
  (net 30 NetR14_1)
  (net 31 NetR12_1)
  (net 32 NetD2_1)
  (net 33 NetD1_1)
  (net 34 NetC5_2)
  (net 35 K_LED4)
  (net 36 K_LED3)
  (net 37 K_LED2)
  (net 38 K_LED1)
  (net 39 GND)
  (net 40 ALS_IN)

  (net_class Default "This is the default net class."
    (add_net ALS_IN)
    (add_net K_LED1)
    (add_net K_LED2)
    (add_net K_LED3)
    (add_net K_LED4)
    (add_net NetC5_2)
    (add_net NetD1_1)
    (add_net NetD2_1)
    (add_net NetR12_1)
    (add_net NetR14_1)
    (add_net NetR15_2)
    (add_net NetR16_2)
    (add_net NetR19_2)
    (add_net NetR3_2)
    (add_net PWM1)
    (add_net PWM2)
    (add_net SK1)
    (add_net SK2)
    (add_net SK3)
    (add_net SK4)
    (add_net S_LED1)
    (add_net S_LED2)
    (add_net S_LED3)
    (add_net S_LED4)
    (add_net S_LED5)
    (add_net S_LED6)
    (add_net S_LED7)
    (add_net TK1)
    (add_net TK2)
    (add_net TK3)
    (add_net TK4)
    (add_net TK5)
    (add_net TK6)
    (add_net VIN)
    (add_net W+)
    (add_net W-)
    (add_net Y+)
    (add_net Y-)
    (trace_width 0.254)
  )
  (net_class Power "This is the default net class."
    (add_net GND)
    (add_net VCC_5V)
    (trace_width 0.8)
  )

 
  
  (module "MIJI_ADPcbLib.PcbLib:R1206-JP_4" (layer Bottom) (tedit 0) (tstamp 0)
    (at 127.600901 96.36789899999977 270)
    (fp_text reference 3 (at 0 0) (layer B.SilkS) hide
      (effects (font (size 0.8 0.8) (thickness 0.127)) (justify left bottom))
    )
    
    (pad 1 smd rect (at 0 0 0) (size 1.5 1) (layers Bottom B.Paste B.Mask)
(net 3 NetR3_2))
     )

  (module "MIJI_ADPcbLib.PcbLib:SMD-EC6.3X6.3X10_5" (layer Bottom) (tedit 0) (tstamp 0)
    (at 114.551001 88.57789899999983 270)
    (fp_text reference 4 (at 0 0) (layer B.SilkS) hide
      (effects (font (size 0.8 0.8) (thickness 0.127)) (justify left bottom))
    )
       (pad 1 smd rect (at 0 0 0) (size 3.3 2.132) (layers Bottom B.Paste B.Mask)
(net 3 NetR3_2))
      )

  
  (module "Miscellaneous Devices LC.PcbLib:0805_C_7" (layer Bottom) (tedit 0) (tstamp 0)
    (at 109.04720099999992 105.11449900000001 270)
    (fp_text reference 6 (at 0 0) (layer B.SilkS) hide
      (effects (font (size 0.8 0.8) (thickness 0.127)) (justify left bottom))
    )
      (pad 1 smd rect (at 0 0 0) (size 1.16 1.47) (layers Bottom B.Paste B.Mask)
(net 6 VIN))

  )

  (module "Miscellaneous Devices LC.PcbLib:1206_C_8" (layer Bottom) (tedit 0) (tstamp 0)
    (at 138.290501 96.60859899999969 180)
    (fp_text reference 7 (at 0 0) (layer B.SilkS) hide
      (effects (font (size 0.8 0.8) (thickness 0.127)) (justify left bottom))
    )
       (pad 1 smd rect (at 0 0 0) (size 1.13 1.8) (layers Bottom B.Paste B.Mask)
(net 6 VIN))
     )

 
  (gr_line (start 160.5011 75.003598) (end 98.501099 75.003598) (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 98.501099 75.003598) (end 98.501099 120.003599) (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 98.501099 120.003599) (end 160.5011 120.003599) (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 160.5011 120.003599) (end 160.5011 75.003598) (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))
))";

    // 两个module的模板
    std::string template_two_modules = R"((kicad_pcb (version 20171130) (host pcbnew "(5.1.2)-2")

  (general
    (thickness 1.6)
  )

  (page A4)
  (layers
    (0 Top signal)
    (31 Bottom signal)
    (32 B.Adhes user)
    (33 F.Adhes user)
    (34 B.Paste user)
    (35 F.Paste user)
    (36 B.SilkS user)
    (37 F.SilkS user)
    (38 B.Mask user)
    (39 F.Mask user)
    (40 Dwgs.User user)
    (41 Cmts.User user)
    (42 Eco1.User user)
    (43 Eco2.User user)
    (44 Edge.Cuts user)
    (45 Margin user)
    (46 B.CrtYd user)
    (47 F.CrtYd user)
    (48 B.Fab user)
    (49 F.Fab user)
  )

  (setup
    (pad_to_mask_clearance 0.1016)
    (aux_axis_origin -120.62617 153.22756)
    (grid_origin -120.62617 153.22756)
    (pcbplotparams
      (layerselection 0x00010fc_ffffffff)
      (plot_on_all_layers_selection 0x0000000_00000000)
      (disableapertmacros false)
      (usegerberextensions false)
      (usegerberattributes true)
      (usegerberadvancedattributes true)
      (creategerberjobfile true)
      (dashed_line_dash_ratio 12)
      (dashed_line_gap_ratio 3)
      (svgprecision 4)
      (plotframeref false)
      (viasonmask false)
      (mode 1)
      (useauxorigin false)
      (hpglpennumber 1)
      (hpglpenspeed 20)
      (hpglpendiameter 15)
      (dxfpolygonmode true)
      (dxfimperialunits true)
      (dxfusepcbnewfont true)
      (psnegative false)
      (psa4output false)
      (plotreference true)
      (plotvalue true)
      (plotinvisibletext false)
      (sketchpadsonfab false)
      (subtractmaskfromsilk false)
      (outputformat 1)
      (mirror false)
      (drillshape 1)
      (scaleselection 1)
      (outputdirectory "")
    )
  )

  (net 0 "")

  (net 1 Y+)
  (net 2 Y-)
  (net 3 NetR3_2)
  (net 4 W-)
  (net 5 W+)
  (net 6 VIN)
  (net 7 VCC_5V)
  (net 8 TK6)
  (net 9 TK5)
  (net 10 TK4)
  (net 11 TK3)
  (net 12 TK2)
  (net 13 TK1)
  (net 14 SK4)
  (net 15 SK3)
  (net 16 SK2)
  (net 17 SK1)
  (net 18 S_LED7)
  (net 19 S_LED6)
  (net 20 S_LED5)
  (net 21 S_LED4)
  (net 22 S_LED3)
  (net 23 S_LED2)
  (net 24 S_LED1)
  (net 25 PWM2)
  (net 26 PWM1)
  (net 27 NetR19_2)
  (net 28 NetR16_2)
  (net 29 NetR15_2)
  (net 30 NetR14_1)
  (net 31 NetR12_1)
  (net 32 NetD2_1)
  (net 33 NetD1_1)
  (net 34 NetC5_2)
  (net 35 K_LED4)
  (net 36 K_LED3)
  (net 37 K_LED2)
  (net 38 K_LED1)
  (net 39 GND)
  (net 40 ALS_IN)

  (net_class Default "This is the default net class."
    (add_net ALS_IN)
    (add_net K_LED1)
    (add_net K_LED2)
    (add_net K_LED3)
    (add_net K_LED4)
    (add_net NetC5_2)
    (add_net NetD1_1)
    (add_net NetD2_1)
    (add_net NetR12_1)
    (add_net NetR14_1)
    (add_net NetR15_2)
    (add_net NetR16_2)
    (add_net NetR19_2)
    (add_net NetR3_2)
    (add_net PWM1)
    (add_net PWM2)
    (add_net SK1)
    (add_net SK2)
    (add_net SK3)
    (add_net SK4)
    (add_net S_LED1)
    (add_net S_LED2)
    (add_net S_LED3)
    (add_net S_LED4)
    (add_net S_LED5)
    (add_net S_LED6)
    (add_net S_LED7)
    (add_net TK1)
    (add_net TK2)
    (add_net TK3)
    (add_net TK4)
    (add_net TK5)
    (add_net TK6)
    (add_net VIN)
    (add_net W+)
    (add_net W-)
    (add_net Y+)
    (add_net Y-)
    (trace_width 0.254)
  )
  (net_class Power "This is the default net class."
    (add_net GND)
    (add_net VCC_5V)
    (trace_width 0.8)
  )


 
  
  (module "MIJI_ADPcbLib.PcbLib:R1206-JP_4" (layer Bottom) (tedit 0) (tstamp 0)
    (at 127.600901 96.36789899999977 270)
    (fp_text reference 3 (at 0 0) (layer B.SilkS) hide
      (effects (font (size 0.8 0.8) (thickness 0.127)) (justify left bottom))
    )
    
    (pad 1 smd rect (at 0 0 0) (size 1.5 1) (layers Bottom B.Paste B.Mask)
(net 3 NetR3_2))
     )

  (module "MIJI_ADPcbLib.PcbLib:SMD-EC6.3X6.3X10_5" (layer Bottom) (tedit 0) (tstamp 0)
    (at 114.551001 88.57789899999983 270)
    (fp_text reference 4 (at 0 0) (layer B.SilkS) hide
      (effects (font (size 0.8 0.8) (thickness 0.127)) (justify left bottom))
    )
       (pad 1 smd rect (at 0 0 0) (size 3.3 2.132) (layers Bottom B.Paste B.Mask)
(net 3 NetR3_2))
      )

 
  (gr_line (start 160.5011 75.003598) (end 98.501099 75.003598) (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 98.501099 75.003598) (end 98.501099 120.003599) (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 98.501099 120.003599) (end 160.5011 120.003599) (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))
  (gr_line (start 160.5011 120.003599) (end 160.5011 75.003598) (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))
))";

    // 处理module坐标替换
    std::vector<PointWithNet> replacementPoints;

    if (modulePoints.size() == 2) {
        // 情况1: 只有2个点，直接使用这两个点替换前两个module
        replacementPoints = modulePoints;
        std::cout << "使用2个点替换前两个module" << std::endl;
    }
    else if (modulePoints.size() >= 4) {
        // 情况2: 有4个或更多点，按网络分组处理
        // 找到有两个相同网络的点
        std::map<std::string, std::vector<PointWithNet>> pointsByNet;
        for (const auto& point : modulePoints) {
            pointsByNet[point.net].push_back(point);
        }

        bool foundMatchingNet = false;
        for (const auto& [net, netPoints] : pointsByNet) {
            if (netPoints.size() >= 2) {
                // 使用前两个相同网络的点
                replacementPoints.push_back(netPoints[0]);
                replacementPoints.push_back(netPoints[1]);
                foundMatchingNet = true;
                std::cout << "找到网络 '" << net << "' 的2个点用于前两个module" << std::endl;
                break;
            }
        }

        // 如果没有找到相同网络的点，使用前4个点
        if (!foundMatchingNet) {
            for (size_t i = 0; i < std::min(modulePoints.size(), size_t(4)); i++) {
                replacementPoints.push_back(modulePoints[i]);
            }
            std::cout << "没有找到相同网络的点，使用前4个点" << std::endl;
        }

        // 添加剩余的点
        size_t addedCount = replacementPoints.size();
        for (const auto& point : modulePoints) {
            if (addedCount >= 4) break;

            // 检查是否已经添加过
            bool alreadyAdded = false;
            for (const auto& addedPoint : replacementPoints) {
                if (point.x == addedPoint.x && point.y == addedPoint.y && point.net == addedPoint.net) {
                    alreadyAdded = true;
                    break;
                }
            }

            if (!alreadyAdded) {
                replacementPoints.push_back(point);
                addedCount++;
            }
        }
    }
    else {
        // 其他情况: 使用所有点
        replacementPoints = modulePoints;
        std::cout << "使用所有 " << modulePoints.size() << " 个点" << std::endl;
    }

    // 选择模板
    std::string template_str;
    if (replacementPoints.size() == 2) {
        template_str = template_two_modules;
        std::cout << "使用两个module的模板" << std::endl;
    }
    else {
        template_str = template_four_modules;
        std::cout << "使用四个module的模板" << std::endl;
    }

    // 处理模板替换
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);

    // 分割模板并按行处理
    std::istringstream templateStream(template_str);
    std::string line;
    int moduleCount = 0;
    int grLineCount = 0;

    while (std::getline(templateStream, line)) {
        // 处理module坐标替换
        if (line.find("(module") != std::string::npos) {
            // 输出module行
            oss << line << "\n";

            // 读取下一行(应该包含at坐标)
            if (std::getline(templateStream, line)) {
                if (line.find("(at") != std::string::npos) {
                    // 替换坐标
                    if (moduleCount < replacementPoints.size()) {
                        const auto& point = replacementPoints[moduleCount];
                        oss << "    (at " << point.x << " " << point.y << " 270)\n";
                        std::cout << "替换module " << (moduleCount + 1) << " 坐标为: ("
                            << point.x << ", " << point.y << "), 网络: " << point.net << std::endl;
                    }
                    else {
                        // 如果没有足够的点，保留原坐标
                        oss << line << "\n";
                    }
                    moduleCount++;
                }
                else {
                    oss << line << "\n";
                }
            }
        }
        // 处理gr_line坐标替换
        else if (line.find("(gr_line") != std::string::npos && expandedCorners.size() == 4) {
            if (grLineCount == 0) {
                // 第一条线: 上边 (LT -> RT)
                oss << "  (gr_line (start " << expandedCorners[0].first << " " << expandedCorners[0].second
                    << ") (end " << expandedCorners[1].first << " " << expandedCorners[1].second
                    << ") (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))\n";
            }
            else if (grLineCount == 1) {
                // 第二条线: 右边 (RT -> RB)
                oss << "  (gr_line (start " << expandedCorners[1].first << " " << expandedCorners[1].second
                    << ") (end " << expandedCorners[2].first << " " << expandedCorners[2].second
                    << ") (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))\n";
            }
            else if (grLineCount == 2) {
                // 第三条线: 下边 (RB -> LB)
                oss << "  (gr_line (start " << expandedCorners[2].first << " " << expandedCorners[2].second
                    << ") (end " << expandedCorners[3].first << " " << expandedCorners[3].second
                    << ") (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))\n";
            }
            else if (grLineCount == 3) {
                // 第四条线: 左边 (LB -> LT)
                oss << "  (gr_line (start " << expandedCorners[3].first << " " << expandedCorners[3].second
                    << ") (end " << expandedCorners[0].first << " " << expandedCorners[0].second
                    << ") (layer Edge.Cuts) (width 0.05) (tstamp 42A38F4))\n";
            }
            grLineCount++;
        }
        else {
            oss << line << "\n";
        }


    }

    // 写入文件
    std::ofstream file(outputFilename);
    if (!file.is_open()) {
        throw std::runtime_error("无法创建输出文件: " + outputFilename);
    }

    file << oss.str();
    file.close();

    std::cout << "完整PCB模板文件已生成: " << outputFilename << std::endl;
    std::cout << "共输出 " << moduleCount << " 个module" << std::endl;
}




// ========== 全局变量存储最新生成的文件名 ==========
std::string g_latestGeneratedFile = "";





//-------------------------------------------------找到矩形外的线段部分-------------------------------------------------//
#include <sstream>

// —— 增强版：带 width 的段信息解析（start/end/layer/net/width）
static bool GetSegmentInfoWithWidth(const std::shared_ptr<Node>& segment,
    double& x1, double& y1, double& x2, double& y2,
    double& width, std::string& layer, std::string& net)
{
    x1 = y1 = x2 = y2 = 0.0; width = 0.0; layer.clear(); net.clear();
    bool hasStart = false, hasEnd = false, hasLayer = false, hasNet = false;

    if (!segment || segment->name != "segment") return false;

    for (const auto& ch : segment->children) {
        if (!ch) continue;
        if (ch->name == "start" && ch->parameters.size() >= 2) {
            x1 = std::stod(ch->parameters[0]);
            y1 = std::stod(ch->parameters[1]);
            hasStart = true;
        }
        else if (ch->name == "end" && ch->parameters.size() >= 2) {
            x2 = std::stod(ch->parameters[0]);
            y2 = std::stod(ch->parameters[1]);
            hasEnd = true;
        }
        else if (ch->name == "layer" && !ch->parameters.empty()) {
            layer = ch->parameters[0];
            hasLayer = true;
        }
        else if (ch->name == "net" && !ch->parameters.empty()) {
            net = ch->parameters[0]; // net id 作为字符串存
            hasNet = true;
        }
        else if (ch->name == "width" && !ch->parameters.empty()) {
            width = std::stod(ch->parameters[0]);
        }
    }
    return hasStart && hasEnd && hasLayer && hasNet;
}

// —— 线参数 t（x = x1 + t*dx, y = y1 + t*dy），选绝对值更大的分量算，避免除零
static inline double ParamT(double x1, double y1, double x2, double y2, double xi, double yi)
{
    const double dx = x2 - x1, dy = y2 - y1;
    if (std::fabs(dx) >= std::fabs(dy)) return (dx == 0.0) ? 0.0 : (xi - x1) / dx;
    return (dy == 0.0) ? 0.0 : (yi - y1) / dy;
}

static inline bool IsZeroLen(double x1, double y1, double x2, double y2, double eps = 1e-9) {
    return std::hypot(x2 - x1, y2 - y1) < eps;
}

// —— 对一个线段与“第一矩形”做裁剪：返回 0~2 段“矩形外部”的子段
//    规则：
//    - 若恰好一端在内一端在外：保留“交点→外端”的那一段（1 段）
//    - 若两端都在外且穿过矩形（通常 2 交点）：保留 [start→t_min] 与 [t_max→end] 两段（2 段）
//    - 若仅擦边/切点（1 交点）：也会得到两段，但零长度段会被过滤
static std::vector<std::array<double, 4>>
ClipToOutsideSegments(const RealRectMM& rect,
    double x1, double y1, double x2, double y2)
{
    std::vector<std::array<double, 4>> outs;
    auto interPts = ComputeLineRectIntersections(x1, y1, x2, y2, rect);
    const bool in1 = IsPointInsideRect(x1, y1, rect);
    const bool in2 = IsPointInsideRect(x2, y2, rect);

    if (interPts.empty()) {
        // 无交点：要么全在外且不相交（忽略），要么全在内（忽略）
        return outs;
    }

    // 计算每个交点的 t，并去重（按 t 去重，防止角点重复）
    struct TI { double t, xi, yi; };
    std::vector<TI> tis;
    tis.reserve(interPts.size());
    for (auto& p : interPts) {
        double t = ParamT(x1, y1, x2, y2, p.first, p.second);
        if (t >= -1e-9 && t <= 1.0 + 1e-9) {
            tis.push_back({ std::clamp(t,0.0,1.0), p.first, p.second });
        }
    }
    if (tis.empty()) return outs;
    std::sort(tis.begin(), tis.end(), [](const TI& a, const TI& b) { return a.t < b.t; });
    // 去重容差
    std::vector<TI> uniq;
    for (const auto& it : tis) {
        if (uniq.empty() || std::fabs(it.t - uniq.back().t) > 1e-9) uniq.push_back(it);
    }

    // —— 情况 A：一端内一端外 → 取距离“内端”最近的交点，保留外部单段
    if (in1 ^ in2) {
        // 找离“内端”最近的交点
        double bestD = 1e100, xi = 0, yi = 0;
        if (in1) {
            for (auto& u : uniq) {
                double d = std::hypot(u.xi - x1, u.yi - y1);
                if (d < bestD) { bestD = d; xi = u.xi; yi = u.yi; }
            }
            // 保留 [交点→终点]
            if (!IsZeroLen(xi, yi, x2, y2)) outs.push_back({ xi,yi,x2,y2 });
        }
        else {
            for (auto& u : uniq) {
                double d = std::hypot(u.xi - x2, u.yi - y2);
                if (d < bestD) { bestD = d; xi = u.xi; yi = u.yi; }
            }
            // 保留 [起点→交点]
            if (!IsZeroLen(x1, y1, xi, yi)) outs.push_back({ x1,y1,xi,yi });
        }
        return outs;
    }

    // —— 情况 B：两端都在外（可能 1 或 2 个交点）
    // 取在 [0,1] 范围内的最小/最大 t
    double tmin = uniq.front().t;
    double tmax = uniq.back().t;

    // 段一： [start(0) → tmin]
    double ax1 = x1, ay1 = y1;
    double ax2 = x1 + (x2 - x1) * tmin, ay2 = y1 + (y2 - y1) * tmin;
    if (!IsZeroLen(ax1, ay1, ax2, ay2)) outs.push_back({ ax1,ay1,ax2,ay2 });

    // 段二： [tmax → end(1)]
    double bx1 = x1 + (x2 - x1) * tmax, by1 = y1 + (y2 - y1) * tmax;
    double bx2 = x2, by2 = y2;
    if (!IsZeroLen(bx1, by1, bx2, by2)) outs.push_back({ bx1,by1,bx2,by2 });

    return outs;
}

// —— 主函数：构建“vector< vector<string> >”，每个子 vector 依次为：
//    [startX, startY, endX, endY, width, layer, net]
//    —— 已支持“直穿矩形”的两段外部线（都会收进来）
std::vector<std::vector<std::string>>
BuildClippedSegmentsOutsideFirstRect(const Grid& grid, int layerId, int blockedAreaId, int expandPixels)
{
    std::vector<std::vector<std::string>> out;

    // 1) 计算“第一矩形”（像素域扩张后的现实外沿）
    auto corners = GetExpandedPCBRealCorners(grid, layerId, blockedAreaId, expandPixels); // LT,RT,RB,LB
    if (corners.size() != 4) return out;
    RealRectMM firstRect{
        corners[0].first, // x_left
        corners[0].second,// y_top
        corners[1].first, // x_right
        corners[2].second // y_bot
    };

    // 2) 只看目标层
    const std::string targetLayer = getLayerNameById(grid, layerId);

    // 3) 扫描并裁剪
    for (const auto& seg : grid.segments) {
        double x1, y1, x2, y2, w; std::string layer, net;
        if (!GetSegmentInfoWithWidth(seg, x1, y1, x2, y2, w, layer, net)) continue;
        if (layer != targetLayer) continue;

        // 与边框相交？（含端点在边上）
        auto inters = ComputeLineRectIntersections(x1, y1, x2, y2, firstRect);
        if (inters.empty()) continue;

        // 裁成 0~2 段“矩形外部”子段
        auto pieces = ClipToOutsideSegments(firstRect, x1, y1, x2, y2);
        if (pieces.empty()) continue;

        // 打包为 vector<string>
        std::ostringstream ss;
        auto pack = [&](double v) { ss.str(""); ss.clear(); ss << std::setprecision(12) << std::fixed << v; return ss.str(); };

        for (auto& seg4 : pieces) {
            std::vector<std::string> rec;
            rec.reserve(7);
            rec.push_back(pack(seg4[0]));
            rec.push_back(pack(seg4[1]));
            rec.push_back(pack(seg4[2]));
            rec.push_back(pack(seg4[3]));
            rec.push_back(pack(w));
            rec.push_back(layer);
            rec.push_back(net);
            out.emplace_back(std::move(rec));
        }
    }
    return out;
}




//-----------------------------------------------------------存储顶角像素---------------------------------------
#include "separate3.h"  // 或者 #include "Separate.h"
#include <stdexcept>
#include <vector>
#include <utility>

// ==== 缓存：最近一次计算出的四角像素顶角（LT, RT, RB, LB） ====
static std::vector<std::pair<int, int>> g_lastFinalCornerTopPx1;
/*
// 把“在边界上的现实坐标”稳定映射到小板内部像素（与 main 里的 toPxInside 一致）
static inline std::pair<int, int>
CornerToPxInside(const Grid& grid, double x_mm, double y_mm, int cornerIdx) {
    // cornerIdx: 0=LT, 1=RT, 2=RB, 3=LB
    const double eps = 1.0 / (grid.inputScale * 1024.0); // 远小于1px
    if (cornerIdx == 1 || cornerIdx == 2) x_mm -= eps;   // RT/RB：向左收一点
    if (cornerIdx == 2 || cornerIdx == 3) y_mm -= eps;   // RB/LB：向上收一点
    return RealMMToPixelPx1(grid, x_mm, y_mm);           // 返回( row , col )，1-based
}

// === 实现 1：把已知的最终现实角点(LT,RT,RB,LB) -> 像素顶角(1-based,row,col) ===
std::vector<std::pair<int, int>>
GetFinalCornerTopPixelsInsidePx1(
    const Grid& grid,
    const std::vector<std::pair<double, double>>& finalCornersLT_RT_RB_LB)
{
    if (finalCornersLT_RT_RB_LB.size() != 4)
        throw std::runtime_error(
            "GetFinalCornerTopPixelsInsidePx1: 角点数量必须为4（顺序：LT,RT,RB,LB）");

    std::vector<std::pair<int, int>> out(4);
    for (int i = 0; i < 4; ++i) {
        out[i] = CornerToPxInside(
            grid,
            finalCornersLT_RT_RB_LB[i].first,
            finalCornersLT_RT_RB_LB[i].second,
            i
        );
    }
    g_lastFinalCornerTopPx1 = out; // 更新缓存
    return out;
}


*/

static inline std::pair<int, int>
CornerToPxInside(const Grid& grid, double x_mm, double y_mm, int cornerIdx) {
    // cornerIdx: 0=LT, 1=RT, 2=RB, 3=LB
    const double eps = 1.0 / (grid.inputScale * 1024.0); // 远小于1px

    // 每个角点都向矩形内部微调，确保落在PCB区域内
    switch (cornerIdx) {
    case 0: // LT(左上)：向右下微调
        x_mm += eps; y_mm += eps;
        break;
    case 1: // RT(右上)：向左下微调  
        x_mm -= eps; y_mm += eps;
        break;
    case 2: // RB(右下)：向左上微调
        x_mm -= eps; y_mm -= eps;
        break;
    case 3: // LB(左下)：向右上微调
        x_mm += eps; y_mm -= eps;
        break;
    }
    return RealMMToPixelPx1(grid, x_mm, y_mm); // 返回( row , col )，1-based
}

// === 实现 1：把已知的最终现实角点(LT,RT,RB,LB) -> 像素顶角(1-based,row,col) ===
std::vector<std::pair<int, int>>
GetFinalCornerTopPixelsInsidePx1(
    const Grid& grid,
    const std::vector<std::pair<double, double>>& finalCornersLT_RT_RB_LB)
{
    if (finalCornersLT_RT_RB_LB.size() != 4)
        throw std::runtime_error(
            "GetFinalCornerTopPixelsInsidePx1: 角点数量必须为4（顺序：LT,RT,RB,LB）");

    std::vector<std::pair<int, int>> out(4);
    for (int i = 0; i < 4; ++i) {
        out[i] = CornerToPxInside(
            grid,
            finalCornersLT_RT_RB_LB[i].first,
            finalCornersLT_RT_RB_LB[i].second,
            i
        );
    }
    g_lastFinalCornerTopPx1 = out; // 更新缓存
    return out;
}

// === 实现 2：自动版（内部自己算最终现实角点，再转像素） ===
std::vector<std::pair<int, int>>
GetFinalCornerTopPixelsInsidePx1_Auto(
    const Grid& grid,
    int layerId, int blockedAreaId,
    int expandPixels, int kWholePixels,
    double boardMinX, double boardMinY,
    double boardMaxX, double boardMaxY)
{
    auto finalCorners = GetExpandedPCBCornersSnapAndClip(
        grid, layerId, blockedAreaId,
        expandPixels, kWholePixels,
        boardMinX, boardMinY, boardMaxX, boardMaxY
    );
    return GetFinalCornerTopPixelsInsidePx1(grid, finalCorners);
}

// === 实现 3：读缓存 ===
const std::vector<std::pair<int, int>>&
GetLastFinalCornerTopPixelsInsidePx1()
{
    return g_lastFinalCornerTopPx1;
}

//==================================================================================================================

// === 在 Separate3.cpp 顶部其它 include 之后补充 ===
#include <random>
#include <array>

// --- 生成 UUID v4（xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx）---
static std::string MakeUUIDv4() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);

    std::array<unsigned char, 16> b{};
    for (auto& v : b) v = static_cast<unsigned char>(dist(gen));
    // version (4)
    b[6] = static_cast<unsigned char>((b[6] & 0x0F) | 0x40);
    // variant (10xxxxxx)
    b[8] = static_cast<unsigned char>((b[8] & 0x3F) | 0x80);

    auto hex2 = [](unsigned char v)->std::string {
        static const char* H = "0123456789abcdef";
        std::string s(2, '0');
        s[0] = H[(v >> 4) & 0xF];
        s[1] = H[v & 0xF];
        return s;
        };

    std::string out;
    out.reserve(36);
    for (int i = 0; i < 16; ++i) {
        out += hex2(b[i]);
        if (i == 3 || i == 5 || i == 7 || i == 9) out += '-';
    }
    return out;
}

// --- 去掉可能带来的引号 ---
static std::string StripQuotes(const std::string& s) {
    if (s.size() >= 2 &&
        ((s.front() == '"' && s.back() == '"') ||
            (s.front() == '\'' && s.back() == '\''))) {
        return s.substr(1, s.size() - 2);
    }
    return s;
}

// --- 规范化层名：Top/Bottom -> F.Cu/B.Cu，其它保持不变 ---
static std::string NormalizeLayerName(std::string L) {
    L = StripQuotes(L);
    if (L == "Top")    return "F.Cu";
    if (L == "Bottom") return "B.Cu";
    return L;
}

// === 把 out & updated 里的线段全部写入到 KiCad 7 文件（通过 KiCadParser 的新增接口）===
// 说明：
// - parser：你当前用于“最新生成的 kicad7.0 文件”的解析器/写入器（已有）
// - out：   外矩形外侧的线段，每个元素为 {startX,startY,endX,endY,width,layer,net}
// - updated：从 KiCad5.0 解析出来的线段（struct Segment 在本文件里已有定义）
// 返回：pair{写入的out条数, 写入的updated条数}
static std::pair<int, int>
AddAllSegmentsIntoNewBoard(KiCadParser& parser,
    const std::vector<std::vector<std::string>>& out,
    const std::vector<Segment>& updated)
{
    int addOutCnt = 0, addUpdCnt = 0;

    // 1) out（字符串版）
    for (size_t i = 0; i < out.size(); ++i) {
        const auto& v = out[i];
        if (v.size() < 7) {
            std::cerr << "[Skip out#" << i << "] 字段不足(需要7个)："
                << v.size() << "\n";
            continue;
        }
        try {
            double sx = std::stod(v[0]);
            double sy = std::stod(v[1]);
            double ex = std::stod(v[2]);
            double ey = std::stod(v[3]);
            double w = std::stod(v[4]);
            std::string layer = NormalizeLayerName(v[5]);
            int net = std::stoi(v[6]);
            std::string ts = MakeUUIDv4();

            // 直接调用 KiCadParser 的新增 segment 接口
            // 注意：layer 与 tstamp 作为 C++ 字符串传入即可（序列化时会自动带引号）
            int newId = parser.addSegmentSimple(sx, sy, ex, ey, w, layer, net, ts);
            if (newId > 0) ++addOutCnt;
        }
        catch (const std::exception& e) {
            std::cerr << "[Skip out#" << i << "] 转换失败: " << e.what() << "\n";
        }
    }

    // 2) updated（结构体版）
    for (size_t i = 0; i < updated.size(); ++i) {
        const auto& s = updated[i];
        try {
            std::string layer = NormalizeLayerName(s.layer);
            std::string ts = MakeUUIDv4();
            int newId = parser.addSegmentSimple(
                s.startX, s.startY,
                s.endX, s.endY,
                s.width,
                layer,
                s.net,
                ts
            );
            if (newId > 0) ++addUpdCnt;
        }
        catch (const std::exception& e) {
            std::cerr << "[Skip updated#" << i << "] 写入失败: " << e.what() << "\n";
        }
    }

    return { addOutCnt, addUpdCnt };
}




//--------------------------------------------------主程序入口-------------------------------------------------
int main() {
    try {
        // === 1) 初始化 Grid ===
        std::string filename = "testcase.kicad_pcb";   // 按需修改
        Grid grid;
        grid.setClearance(0.7);
        grid.setInputScale(10);
        grid.SetUp(filename);
        std::cout << "Grid初始化成功! 网格尺寸: " << grid.width << " x " << grid.height << std::endl;

        std::cout << std::fixed << std::setprecision(6) << grid.min_x << ' ' << grid.min_y << ' ' << grid.max_x << ' ' << grid.max_y;

        // === 2) 参数 ===
        int layerId = 1;   // 层序号
        int blockedAreaId = 2;   // 阻塞区域编号
        int expandPixels = 15; //15   // 像素域扩展像素数（有补偿）
        int kWholePixels = 20; //20  // 现实域再扩 k 个整像素（与 0.5px 一起）
        double pxSizeMM = 1.0 / static_cast<double>(grid.inputScale);

        std::cout << "参数: 层=" << layerId << " (" << getLayerNameById(grid, layerId)
            << "), 阻塞区域=" << blockedAreaId
            << ", 像素域扩展=" << expandPixels << " px"
            << ", 现实域最终再扩 = 0.5px + " << kWholePixels << "px"
            << " ≈ " << std::fixed << std::setprecision(5)
            << (0.5 + kWholePixels) * pxSizeMM << " mm\n";

        // === 3) 读取真实板框（Edge.Cuts 外接矩形） ===
        KiCadParser parser2;
        if (!parser2.parseFile(filename)) {
            throw std::runtime_error("解析 KiCad 文件失败，无法获取板框");
        }
        KiCadParser::Point2D bl, tl, tr, br; // 左下、左上、右上、右下
        if (!parser2.getBoardCorners(bl, tl, tr, br)) {
            throw std::runtime_error("未找到 Edge.Cuts 板框");
        }
        const double boardMinX = bl.x;
        const double boardMinY = bl.y;
        const double boardMaxX = tr.x;
        const double boardMaxY = tr.y;

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[Board] min=(" << boardMinX << "," << boardMinY
            << "), max=(" << boardMaxX << "," << boardMaxY << ")\n";

        // === 4) 打印"像素扩展阶段"的像素坐标（有补偿） ===
        // 4.1 阻塞区域像素（1-based）
        auto blockedPixels = GetBlockedAreaPixels_1Based(grid, layerId, blockedAreaId);
        if (blockedPixels.empty()) {
            throw std::runtime_error("未找到指定的阻塞区域");
        }

        // 4.2 原始像素包围盒（1-based）
        BBoxPx1 originalBBox = ComputeBBoxFromPixels1(blockedPixels);
        std::cout << "\n=== 原始像素包围盒(1-based, col,x  row,y) ===\n";
        std::cout << "LT(col,row)=(" << originalBBox.col_min << "," << originalBBox.row_min << ")\n";
        std::cout << "RT(col,row)=(" << originalBBox.col_max << "," << originalBBox.row_min << ")\n";
        std::cout << "RB(col,row)=(" << originalBBox.col_max << "," << originalBBox.row_max << ")\n";
        std::cout << "LB(col,row)=(" << originalBBox.col_min << "," << originalBBox.row_max << ")\n";
        std::cout << "W x H (px) = " << (originalBBox.col_max - originalBBox.col_min + 1)
            << " x " << (originalBBox.row_max - originalBBox.row_min + 1) << "\n";

        // 4.3 像素扩展（有补偿）
        ClipReport clipReport{};
        BBoxPx1 expandedBBoxPx = ExpandBBoxByMPixelsCompensated(originalBBox, expandPixels, grid, &clipReport);

        std::cout << "\n=== 像素扩展后包围盒(有补偿, 1-based) ===\n";
        std::cout << "LT(col,row)=(" << expandedBBoxPx.col_min << "," << expandedBBoxPx.row_min << ")\n";
        std::cout << "RT(col,row)=(" << expandedBBoxPx.col_max << "," << expandedBBoxPx.row_min << ")\n";
        std::cout << "RB(col,row)=(" << expandedBBoxPx.col_max << "," << expandedBBoxPx.row_max << ")\n";
        std::cout << "LB(col,row)=(" << expandedBBoxPx.col_min << "," << expandedBBoxPx.row_max << ")\n";
        std::cout << "W x H (px) = " << (expandedBBoxPx.col_max - expandedBBoxPx.col_min + 1)
            << " x " << (expandedBBoxPx.row_max - expandedBBoxPx.row_min + 1) << "\n";
        std::cout << "触边: left=" << clipReport.left
            << ", right=" << clipReport.right
            << ", top=" << clipReport.top
            << ", bottom=" << clipReport.bottom << "\n";

        // === 5) （可选）打印"像素域扩展后的现实四角"（用于对比）
        auto baseCorners = GetExpandedPCBRealCorners(grid, layerId, blockedAreaId, expandPixels);
        std::vector<std::string> names = { "左上(LT)","右上(RT)","右下(RB)","左下(LB)" };
        std::cout << "\n=== 像素域扩展后的现实四角(用于对比) ===\n";
        for (size_t i = 0; i < baseCorners.size(); ++i) {
            std::cout << names[i] << ": (" << baseCorners[i].first << ", " << baseCorners[i].second << ")\n";
        }

        // === 6) 最终：现实域 0.5px + k×1px 对齐扩张 + 裁剪到板框 ===
        auto expandedCorners = GetExpandedPCBCornersSnapAndClip(
            grid, layerId, blockedAreaId,
            expandPixels, kWholePixels,
            boardMinX, boardMinY, boardMaxX, boardMaxY
        );

        std::cout << "\n=== 最终小PCB边界框(0.5px + " << kWholePixels
            << "px 对齐 + 板框裁剪后) ===\n";

        for (size_t i = 0; i < expandedCorners.size(); ++i) {
            std::cout << names[i] << ": (" << expandedCorners[i].first
                << ", " << expandedCorners[i].second << ")\n";
        }

        // === 6.1) 最终四角 -> 像素顶角坐标(1-based, row,col)，保证落在小板内部 ===
        auto toPxInside = [&](double x_mm, double y_mm, int cornerIdx) -> std::pair<int, int> {
            // cornerIdx: 0=LT, 1=RT, 2=RB, 3=LB
            const double eps = 1.0 / (grid.inputScale * 1024.0); // 远小于1px
            double xx = x_mm, yy = y_mm;

            // 每个角点都向矩形内部微调，确保落在PCB区域内
            switch (cornerIdx) {
            case 0: // LT(左上)：向右下微调
                xx += eps; yy += eps;
                break;
            case 1: // RT(右上)：向左下微调  
                xx -= eps; yy += eps;
                break;
            case 2: // RB(右下)：向左上微调
                xx -= eps; yy -= eps;
                break;
            case 3: // LB(左下)：向右上微调
                xx += eps; yy -= eps;
                break;
            }
            return RealMMToPixelPx1(grid, xx, yy); // 返回( row , col )，1-based
            };

        std::array<std::pair<int, int>, 4> finalCornerPx;
        for (int i = 0; i < 4; ++i) {
            finalCornerPx[i] = toPxInside(expandedCorners[i].first, expandedCorners[i].second, i);
        }

        std::cout << "\n=== 最终小PCB四角对应的像素顶角坐标(1-based, row,col) ===\n";
        const char* nm[4] = { "LT","RT","RB","LB" };
        for (int i = 0; i < 4; ++i) {
            std::cout << nm[i] << ": (" << finalCornerPx[i].first
                << ", " << finalCornerPx[i].second << ")\n";
        }

        // === 6.2) 计算并打印边界框的尺寸 ===
        std::cout << "\n=== 小PCB边界框尺寸计算 ===\n";

        // 物理尺寸（毫米）
        double physicalWidth = expandedCorners[1].first - expandedCorners[0].first;  // RT.x - LT.x
        double physicalHeight = expandedCorners[3].second - expandedCorners[0].second; // LB.y - LT.y

        std::cout << "物理宽度 (RT.x - LT.x): " << expandedCorners[1].first << " - "
            << expandedCorners[0].first << " = " << physicalWidth << " mm\n";
        std::cout << "物理高度 (LB.y - LT.y): " << expandedCorners[3].second << " - "
            << expandedCorners[0].second << " = " << physicalHeight << " mm\n";

        // 像素尺寸
        int pixelWidth = finalCornerPx[1].second - finalCornerPx[0].second;  // RT.col - LT.col
        int pixelHeight = finalCornerPx[3].first - finalCornerPx[0].first;   // LB.row - LT.row

        std::cout << "像素宽度 (RT.col - LT.col): " << finalCornerPx[1].second << " - "
            << finalCornerPx[0].second << " = " << pixelWidth << " px\n";
        std::cout << "像素高度 (LB.row - LT.row): " << finalCornerPx[3].first << " - "
            << finalCornerPx[0].first << " = " << pixelHeight << " px\n";

        // 也可以计算对角线的尺寸
        double physicalWidth2 = expandedCorners[2].first - expandedCorners[3].first;  // RB.x - LB.x
        double physicalHeight2 = expandedCorners[2].second - expandedCorners[1].second; // RB.y - RT.y

        std::cout << "\n验证尺寸（另一组边）:\n";
        std::cout << "物理宽度 (RB.x - LB.x): " << expandedCorners[2].first << " - "
            << expandedCorners[3].first << " = " << physicalWidth2 << " mm\n";
        std::cout << "物理高度 (RB.y - RT.y): " << expandedCorners[2].second << " - "
            << expandedCorners[1].second << " = " << physicalHeight2 << " mm\n";

        int pixelWidth2 = finalCornerPx[2].second - finalCornerPx[3].second;  // RB.col - LB.col
        int pixelHeight2 = finalCornerPx[2].first - finalCornerPx[1].first;   // RB.row - RT.row

        std::cout << "像素宽度 (RB.col - LB.col): " << finalCornerPx[2].second << " - "
            << finalCornerPx[3].second << " = " << pixelWidth2 << " px\n";
        std::cout << "像素高度 (RB.row - RT.row): " << finalCornerPx[2].first << " - "
            << finalCornerPx[1].first << " = " << pixelHeight2 << " px\n";

        // === 使用通用 get 函数（基于已得到的 expandedCorners） ===
        auto finalCornerPxVec = GetFinalCornerTopPixelsInsidePx1(grid, expandedCorners);

        std::cout << "\n=== 最终小PCB四角对应的像素顶角坐标(1-based, row,col) [via GetFinalCornerTopPixelsInsidePx1] ===\n";
        const char* nm2[4] = { "LT","RT","RB","LB" };
        for (int i = 0; i < 4; ++i) {
            std::cout << nm2[i] << ": (" << finalCornerPxVec[i].first
                << ", " << finalCornerPxVec[i].second << ")\n";
        }






        // === 计算像素面积的两倍 ===
        std::cout << "\n=== 小PCB边界框像素面积×2 ===\n";

        // 直接使用之前计算的像素尺寸（不需要重新定义）
        int pixelAreaTimes2 = pixelWidth * pixelHeight * 2;

        std::cout << "像素宽度: " << pixelWidth << " px\n";
        std::cout << "像素高度: " << pixelHeight << " px\n";
        std::cout << "像素面积×2: " << pixelWidth << " × " << pixelHeight << " × 2 = " << pixelAreaTimes2 << " px²\n";







        // === 7) CostUpdater 部分 ===
        CostUpdater costUpdater(grid);
        std::string outputDir = "D:\\梁彦诗的项目\\Liang的科研\\第五步最后9999999999999999999999999999999999999\\output\\";

        // 创建输出目录（如果不存在）
        // 注意：这里需要检查目录是否存在并创建，但为了避免引入额外依赖，我们先假设目录存在
        std::string weightsOutputFile = outputDir + "weights_output.txt";

        costUpdater.updateCost(layerId, finalCornerPxVec[0].second, finalCornerPxVec[0].first,
            finalCornerPxVec[2].second, finalCornerPxVec[2].first);
        costUpdater.outputWeightsToFile(layerId, finalCornerPxVec[0].second, finalCornerPxVec[0].first,
            finalCornerPxVec[2].second, finalCornerPxVec[2].first,
            weightsOutputFile);

        std::cout << "权重文件已输出到: " << weightsOutputFile << std::endl;



        // === 8) 交点/终止点（仍按像素域扩展的矩形做，保持你原逻辑） ===
        auto points = FindIntersectionsAndEndpointsInSmallPCB(grid, layerId, blockedAreaId, expandPixels);


        /*
        std::cout << "\n=== 线段交点和终止点 ===\n";
        std::cout << "共找到 " << points.size() << " 个点\n";

        */


        // === 9) 生成模板（用最终四角）===
        std::string templateOutputFile = outputDir + "complete_pcb_template.kicad_pcb";
        GenerateCompletePCBTemplate(points, expandedCorners, templateOutputFile);
        std::cout << "PCB模板已输出到: " << templateOutputFile << std::endl;

        // === 10) 测试删除的segment并生成新的kicad7.0文件 ===
        // 收集并删除"第一步像素扩张形成的现实外沿矩形"内/边上的线段
        {
            // 1) 收集候选 SEID（注意：expandPixels 必须与第一步像素扩展保持一致）
            std::vector<int> candIds = CollectSegmentsTouchingFirstRect(
                grid,           // 当前 Grid（含 segments / layerName 等）
                layerId,        // 这层
                blockedAreaId,  // 这块阻塞区域
                expandPixels    // 第一步像素扩张用的像素数
            );
            std::sort(candIds.begin(), candIds.end());
            candIds.erase(std::unique(candIds.begin(), candIds.end()), candIds.end());

            std::cout << "\n=== 待删除线段（第一矩形内/边上）===\n";
            if (candIds.empty()) {
                std::cout << "无\n";
            }
            else {
                std::cout << "共 " << candIds.size() << " 条：";
                for (size_t i = 0; i < candIds.size(); ++i) {
                    if (i) std::cout << ", ";
                    std::cout << candIds[i];
                }
                std::cout << "\n";

                // —— 明细打印：SEID + 起终点 + 层 + net
                auto extractSeid = [](const std::shared_ptr<Node>& seg) -> int {
                    if (!seg || seg->parameters.empty()) return -1;
                    const std::string& p0 = seg->parameters.front(); // 形如 [[SEID:123]]
                    const std::string prefix = "[[SEID:";
                    const std::string suffix = "]]";
                    if (p0.rfind(prefix, 0) != 0) return -1;
                    size_t end = p0.rfind(suffix);
                    if (end == std::string::npos) return -1;
                    try { return std::stoi(p0.substr(prefix.size(), end - prefix.size())); }
                    catch (...) { return -1; }
                    };

                std::cout << "—— 明细 ——\n";
                for (const auto& seg : grid.segments) {
                    int sid = extractSeid(seg);
                    if (sid <= 0) continue;
                    if (!std::binary_search(candIds.begin(), candIds.end(), sid)) continue;

                    double x1, y1, x2, y2; std::string lyr, net;
                    if (GetSegmentInfo(seg, x1, y1, x2, y2, lyr, net)) {
                        std::cout << "  SEID " << sid << " [" << lyr << "]: "
                            << "(" << x1 << ", " << y1 << ") -> ("
                            << x2 << ", " << y2 << ")  net=" << net << "\n";
                    }
                }
            }

            // 2) 真正删除并另存新文件
            if (!candIds.empty()) {
                KiCadParser parser;
                if (!parser.parseFile(filename)) {
                    std::cerr << "解析原始 Kicad 文件失败: " << filename << "\n";
                }
                else {
                    size_t removed = 0;
                    for (int id : candIds) removed += parser.removeSegmentById(id);
                    std::cout << "已删除线段条数: " << removed << "\n";

                    // 输出文件名到指定目录
                    std::string outPath = outputDir + "firstrect_pruned.kicad_pcb";

                    if (parser.saveAsKicadPcb(outPath, /*stripTempIds=*/true)) {
                        std::cout << "已写出删除后的文件: " << outPath << "\n";
                    }
                    else {
                        std::cerr << "写出失败\n";
                    }

                    // ====== (新增) 生成"矩形外"的片段集合（vector<vector<string>>）======
                    auto out_vec = BuildClippedSegmentsOutsideFirstRect(
                        grid,        // 7.0 的 grid
                        layerId,
                        blockedAreaId,
                        expandPixels // 与删除步骤一致
                    );
                    std::cout << "[新增-矩形外] 条数: " << out_vec.size() << "\n";

                    // ====== (新增) 从 KiCad 5.0 文件解析得到小板段（seg5）======
                    std::string v5InputPath =
                        "D:\\梁彦诗的项目\\Liang的科研\\第五步最后9999999999999999999999999999999999999\\output\\test.kicad_pcb";

                    std::vector<Segment> seg5;
                    std::vector<SegmentSpan> sp5;
                    std::string v5Raw; // 可选：保留原文本
                    LoadSegmentsWithSpans(v5InputPath, seg5, sp5, &v5Raw);
                    std::cout << "[v5解析] 段数: " << seg5.size() << "\n";

                    // ====== (新增) 从 v7 原始板（filename）解析大板段（seg7），用于匹配net/width ======
                    std::vector<Segment> seg7;
                    std::vector<SegmentSpan> sp7;
                    LoadSegmentsWithSpans(filename, seg7, sp7, nullptr);
                    std::cout << "[v7解析] 段数: " << seg7.size() << "\n";

                    // ====== (新增) 从小板v5文件加载pads ======
                    std::vector<Pad2> pads;
                    LoadPadsAbsolute(v5InputPath, pads);

                    // ====== (新增) 逐段匹配 + 路径投票，修正 seg5 的 net/width ======
                    MatchParams mp;
                    mp.coordTol = 0.05;           // 容差（可调）
                    mp.requireLayerEqual = true;  // 要求同层匹配

                    std::vector<int>                 newNet(seg5.size(), INT_MIN);
                    std::vector<double>              newW(seg5.size(), -1.0);
                    std::vector<std::string>         newLayer(seg5.size(), "");

                    std::vector<std::vector<std::pair<int, double>>> hits(seg5.size());

                    // 先用 pad 做种子（这一步是新增的关键）
                    const double PAD_TOL = 0.25, ATTACH_TOL = 0.25;  // 和 222 保持一致
                    SeedSmallFromPads(pads, seg7 /*原板段*/, seg5 /*小板段*/, PAD_TOL, ATTACH_TOL,
                        newNet, newW, newLayer);

                    // 然后对"还没被 pad 赋值"的小段，走你现有的端点投票：
                    for (int i = 0; i < (int)seg5.size(); ++i) {
                        if (newNet[i] != INT_MIN && newW[i] >= 0) continue;  // 已经有 pad 结果就跳过
                        int net; double w;
                        pickOrigNetAndWidthForSmallSegment(seg5[i], seg7, MatchParams{ 0.05,true }, net, w, nullptr);
                        if (net != -1) { newNet[i] = net; newW[i] = w; newLayer[i] = seg5[i].layer; }
                    }

                    // 逐段匹配（只收集hits，不覆盖newNet/newW）
                    for (size_t i = 0; i < seg5.size(); ++i) {
                        int tnet; double tw;
                        pickOrigNetAndWidthForSmallSegment(seg5[i], seg7, mp, tnet, tw, &hits[i]);
                    }


#if 0
                    // 路径聚类（按端点共享 + 层内）
                    std::vector<std::vector<int>> adj;
                    buildAdjacencyBySharedEndpoints(seg5, mp.coordTol, adj);
                    auto comps = connectedComponents(adj);

                    // 在每条路径内做多数投票，统一该路径所有段的 net/width
                    for (auto& comp : comps) {
                        // 投 net
                        std::map<int, int> voteN;
                        for (int si : comp)
                            for (auto& h : hits[si]) voteN[h.first]++;
                        if (voteN.empty()) continue;

                        int bestNet = -1, bestVotes = -1;
                        for (auto& kv : voteN) {
                            if (kv.second > bestVotes) { bestVotes = kv.second; bestNet = kv.first; }
                        }

                        // 在确定的 net 下投 width
                        std::map<double, int> voteW;
                        for (int si : comp)
                            for (auto& h : hits[si])
                                if (h.first == bestNet) voteW[h.second]++;
                        if (!voteW.empty()) {
                            int wv = -1; double bestW = 0;
                            for (auto& kv : voteW) {
                                if (kv.second > wv) { wv = kv.second; bestW = kv.first; }
                            }
                            for (int si : comp) {
                                newNet[si] = bestNet;
                                newW[si] = bestW;
                            }
                        }
                    }
#endif





                    // —— 目标层名（统一强制到这个层）
                    const std::string targetLayer = getLayerNameById(grid, layerId);

                    // —— 对仍未被 pad/端点命中的段，用“目标层上的最近原段”补 net/width（非投票）
                    for (size_t i = 0; i < seg5.size(); ++i) {
                        if (newNet[i] != INT_MIN && newW[i] >= 0) continue; // 已经有结果的不动

                        // 用段中点做最近匹配（也可以换成端点最小距离）
                        double mx = 0.5 * (seg5[i].startX + seg5[i].endX);
                        double my = 0.5 * (seg5[i].startY + seg5[i].endY);

                        double best = 1e100; int pick = -1;
                        for (int j = 0; j < (int)seg7.size(); ++j) {
                            if (seg7[j].layer != targetLayer) continue;
                            bool on = false;
                            double d2 = pointSegmentDist2(mx, my, seg7[j].startX, seg7[j].startY, seg7[j].endX, seg7[j].endY, on);
                            // 不要求“在线上”，纯最近；如果你想更严，可以加 on 条件
                            if (d2 < best) { best = d2; pick = j; }
                        }
                        if (pick >= 0) {
                            if (newNet[i] == INT_MIN) newNet[i] = seg7[pick].net;
                            if (newW[i] < 0)       newW[i] = seg7[pick].width;
                            // 层先不动，这里只是确定 net/width
                        }
                    }







                    // 形成"修正后的 v5 段集合"，后续将把它们写回 v7
                    /*std::vector<Segment> updated_fixed;
                    updated_fixed.reserve(seg5.size());
                    for (size_t i = 0; i < seg5.size(); ++i) {
                        Segment s = seg5[i];
                        if (newNet[i] != INT_MIN) s.net = newNet[i];
                        if (newW[i] >= 0) s.width = newW[i];
                        if (!newLayer[i].empty()) s.layer = newLayer[i];
                        updated_fixed.push_back(std::move(s));
                    }*/


                    // 形成"修正后的 v5 段集合"
                    std::vector<Segment> updated_fixed;
                    updated_fixed.reserve(seg5.size());
                    for (size_t i = 0; i < seg5.size(); ++i) {
                        Segment s = seg5[i];
                        if (newNet[i] != INT_MIN) s.net = newNet[i];
                        if (newW[i] >= 0)        s.width = newW[i];

                        // ✅ 强制所有小段写入指定层（按你传入的 layerId）
                        s.layer = targetLayer;  // 等同于 getLayerNameById(grid, layerId)

                        updated_fixed.push_back(std::move(s));
                    }

                    // （可选）打印 v5(修正后) 将要加入的线段（前10条）
                    std::cout << "[V5(修正后) 需新增线段] 共 " << updated_fixed.size() << " 条\n";
                    for (size_t i = 0; i < std::min<size_t>(10, updated_fixed.size()); ++i) {
                        const auto& s = updated_fixed[i];
                        std::cout << "  #" << i
                            << " start(" << s.startX << ", " << s.startY << ")"
                            << " end(" << s.endX << ", " << s.endY << ")"
                            << " width=" << s.width
                            << " layer=" << s.layer
                            << " net=" << s.net
                            << "\n";
                    }

                    // ====== (新增) 读入"删除后的文件"，把两类 segment 写进去（矩形外 + v5修正后）======
                    KiCadParser parserAdd;
                    if (!parserAdd.parseFile(outPath)) { // outPath 是你上一步保存的 ".firstrect_pruned.kicad_pcb"
                        std::cerr << "无法读取删除后的文件: " << outPath << "\n";
                        return 2;
                    }

                    auto addStat = AddAllSegmentsIntoNewBoard(parserAdd, out_vec, updated_fixed);
                    std::cout << "[新增执行] 外矩形外加入: " << addStat.first
                        << " 条；v5(已修正)加入: " << addStat.second << " 条\n";

                    // ====== (新增) 再次保存 —— 终版 v7 文件（存到指定输出目录）======
                    std::string outPath2 = outputDir + "merged_final_v7.kicad_pcb";

                    if (parserAdd.saveAsKicadPcb(outPath2, /*stripTempIds=*/true)) {
                        std::cout << "已写出最终文件: " << outPath2 << "\n";
                    }
                    else {
                        std::cerr << "最终文件写出失败\n";
                    }
                }
            }
        }

        // === [TEST v2] 原始 segment vs. 裁剪后的外部片段（含直穿矩形的两段） ===
        {
            // 计算第一矩形（像素扩张后的现实外沿）——若你外面已有 corners，可直接复用
            auto corners = GetExpandedPCBRealCorners(grid, layerId, blockedAreaId, expandPixels); // LT, RT, RB, LB
            if (corners.size() != 4) {
                std::cerr << "[TEST v2] 第一矩形计算失败，corners.size()=" << corners.size() << "\n";
            }
            else {
                RealRectMM firstRect{
                    corners[0].first, // x_left
                    corners[0].second,// y_top
                    corners[1].first, // x_right
                    corners[2].second // y_bot
                };

                const std::string targetLayer = getLayerNameById(grid, layerId);
                size_t totalPieces = 0, touchedSegs = 0;

                std::cout << "\n=== 原始 segment 与裁剪后片段对应关系（层: " << targetLayer << "） ===\n";
                std::cout << std::fixed << std::setprecision(6);

                for (const auto& seg : grid.segments) {
                    double x1, y1, x2, y2, w; std::string layer, net;
                    if (!GetSegmentInfoWithWidth(seg, x1, y1, x2, y2, w, layer, net)) continue;
                    if (layer != targetLayer) continue;

                    // 只看与矩形边框有交点的 segment（端点在边上也算）
                    auto inters = ComputeLineRectIntersections(x1, y1, x2, y2, firstRect);
                    if (inters.empty()) continue;

                    // 把"矩形内部"部分剔除，保留矩形外的 0~2 段
                    auto pieces = ClipToOutsideSegments(firstRect, x1, y1, x2, y2);
                    if (pieces.empty()) continue;

                    ++touchedSegs;
                    std::cout << "原始: (" << x1 << ", " << y1 << ") -> ("
                        << x2 << ", " << y2 << ")  w=" << w
                        << ", layer=" << layer << ", net=" << net << "\n";

                    for (size_t k = 0; k < pieces.size(); ++k) {
                        const auto& p = pieces[k]; // [sx, sy, ex, ey]
                        std::cout << "    裁剪[" << k << "]: ("
                            << p[0] << ", " << p[1] << ") -> ("
                            << p[2] << ", " << p[3] << ")\n";
                        ++totalPieces;
                    }
                }

                // 顺便跑一下"打包后的新 vector"，确认数量一致
                auto packed = BuildClippedSegmentsOutsideFirstRect(
                    grid,        // Grid（含 segments）
                    layerId,     // 目标层
                    blockedAreaId,
                    expandPixels // 必须与第一步像素扩张一致
                );

                std::cout << "\n=== 汇总 ===\n";
                std::cout << "涉及到的原始 segment 条数: " << touchedSegs << "\n";
                std::cout << "裁剪得到的外部片段总数:   " << totalPieces << "\n";
                std::cout << "打包 vector 的条目数:     " << packed.size() << "\n";


                // === 计算像素面积的两倍 ===
                std::cout << "\n=== 小PCB边界框像素面积×2 ===\n";

                // 直接使用之前计算的像素尺寸（不需要重新定义）
                int pixelAreaTimes2 = pixelWidth * pixelHeight * 2;

                std::cout << "像素宽度: " << pixelWidth << " px\n";
                std::cout << "像素高度: " << pixelHeight << " px\n";
                std::cout << "像素面积×2: " << pixelWidth << " × " << pixelHeight << " × 2 = " << pixelAreaTimes2 << " px²\n";

            }
        }

    }
    catch (const std::exception& e) {
        std::cerr << "程序运行出错: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "程序执行完成!" << std::endl;
    return 0;
}