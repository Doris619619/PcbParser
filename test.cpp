/*
PCB连通性图构建算法
充分利用Lemon图库内置算法进行全局图分析
编译: g++ test.cpp -o test -std=c++17 -I$HOME/local/include -L$HOME/local/lib -lemon
*/

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <memory>
#include <queue>
#include <fstream>
#include <climits>   

// Lemon库头文件
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/bfs.h>
#include <lemon/dijkstra.h>
#include <lemon/dfs.h>
#include <lemon/path.h>
#include <lemon/kruskal.h>

using namespace std;
using namespace lemon;

// 区域结构定义 
struct Region {
    int id;
    string type;                               // GND, VCC, Signal, Blocked
    int layer;
    vector<pair<int, int>> pixels;
    
    struct BoundingBox {
        int minX, maxX, minY, maxY;
        BoundingBox() : minX(INT_MAX), maxX(INT_MIN), minY(INT_MAX), maxY(INT_MIN) {}
    } boundingBox;
    
    Region(int _id, const string& _type, int _layer) 
        : id(_id), type(_type), layer(_layer) {}
    
    void addPixel(int x, int y) {
        pixels.push_back({x, y});
        boundingBox.minX = min(boundingBox.minX, x);
        boundingBox.maxX = max(boundingBox.maxX, x);
        boundingBox.minY = min(boundingBox.minY, y);
        boundingBox.maxY = max(boundingBox.maxY, y);
    }
    
    pair<int, int> getCenter() const {
        if (pixels.empty()) return {0, 0};
        int sumX = 0, sumY = 0;
        for (const auto& pixel : pixels) {
            sumX += pixel.first;
            sumY += pixel.second;
        }
        return {sumX / pixels.size(), sumY / pixels.size()};
    }
};

struct Cell {
    int layer;
    int x, y;
    bool isBlocked;
    bool isEmpty;
    shared_ptr<Region> regionPtr;
    
    Cell() : layer(0), x(0), y(0), isBlocked(false), isEmpty(true), regionPtr(nullptr) {}
    Cell(int _layer, int _x, int _y) : layer(_layer), x(_x), y(_y), isBlocked(false), isEmpty(true), regionPtr(nullptr) {}
};

struct PathFindResult {
    bool canConnect;
    int blockedPixelsToRemove;
    vector<pair<int, int>> pathPixels;
    string connectionType;
    double pathCost;
    
    PathFindResult() : canConnect(false), blockedPixelsToRemove(0), connectionType("impossible"), pathCost(0.0) {}
};

//PCB连通性图构建器
class PCBConnectivityGraphBuilder {
private:
    ListGraph graph;
    unique_ptr<ListGraph::NodeMap<int>> nodeRegionId;
    unique_ptr<ListGraph::NodeMap<string>> nodeType;     
    unique_ptr<ListGraph::NodeMap<int>> nodeLayer;       
    unique_ptr<ListGraph::EdgeMap<int>> edgeWeight;
    unique_ptr<ListGraph::EdgeMap<string>> edgeType;
    unique_ptr<ListGraph::EdgeMap<double>> edgeCost;     
    
    map<int, ListGraph::Node> regionToNode;
    vector<Region> regions;
    vector<vector<vector<Cell>>> pcbMatrix;
    
    int overlapThreshold;
    int maxClearableBlocked;

public:
    PCBConnectivityGraphBuilder(int threshold = 3, int maxBlocked = 50) 
        : overlapThreshold(threshold), maxClearableBlocked(maxBlocked) {
        initializeMaps();
    }
    
    void initializeMaps() {
        // 使用智能指针管理Map对象，避免赋值问题
        nodeRegionId = make_unique<ListGraph::NodeMap<int>>(graph);
        nodeType = make_unique<ListGraph::NodeMap<string>>(graph);
        nodeLayer = make_unique<ListGraph::NodeMap<int>>(graph);
        edgeWeight = make_unique<ListGraph::EdgeMap<int>>(graph);
        edgeType = make_unique<ListGraph::EdgeMap<string>>(graph);
        edgeCost = make_unique<ListGraph::EdgeMap<double>>(graph);
    }
    
    void extractRegionsFromMatrix(const vector<vector<vector<Cell>>>& matrix) {
        pcbMatrix = matrix;
        regions.clear();
        map<string, shared_ptr<Region>> regionMap;
        int regionIdCounter = 0;
        
        for (int layer = 0; layer < matrix.size(); layer++) {
            for (int x = 0; x < matrix[layer].size(); x++) {
                for (int y = 0; y < matrix[layer][x].size(); y++) {
                    const Cell& cell = matrix[layer][x][y];
                    
                    if (!cell.isEmpty && !cell.isBlocked && cell.regionPtr) {
                        string regionKey = cell.regionPtr->type + "_L" + to_string(layer) + 
                                         "_" + to_string(x/5) + "_" + to_string(y/5);
                        
                        if (regionMap.find(regionKey) == regionMap.end()) {
                            regionMap[regionKey] = make_shared<Region>(regionIdCounter++, cell.regionPtr->type, layer);
                        }
                        
                        regionMap[regionKey]->addPixel(x, y);
                    }
                }
            }
        }
        
        for (auto& pair : regionMap) {
            regions.push_back(*pair.second);
        }
        
        cout << "提取到 " << regions.size() << " 个区域" << endl;
    }
    
    void buildConnectivityGraph() {
        // 重新构建图和所有映射
        graph.clear();
        regionToNode.clear();
        initializeMaps();
        
        // 创建节点
        for (const auto& region : regions) {
            ListGraph::Node node = graph.addNode();
            (*nodeRegionId)[node] = region.id;
            (*nodeType)[node] = region.type;
            (*nodeLayer)[node] = region.layer;
            regionToNode[region.id] = node;
        }
        
        cout << "创建了 " << countNodes(graph) << " 个节点" << endl;
        
        // 创建边
        int directEdges = 0, blockedEdges = 0, crossLayerEdges = 0;
        
        for (int i = 0; i < regions.size(); i++) {
            for (int j = i + 1; j < regions.size(); j++) {
                PathFindResult result = analyzeRegionConnection(regions[i], regions[j]);
                
                if (result.canConnect) {
                    ListGraph::Node u = regionToNode[regions[i].id];
                    ListGraph::Node v = regionToNode[regions[j].id];
                    
                    if (findEdge(graph, u, v) == INVALID) {
                        ListGraph::Edge edge = graph.addEdge(u, v);
                        
                        if (result.connectionType == "direct") {
                            (*edgeWeight)[edge] = 1;
                            (*edgeCost)[edge] = 1.0;
                            (*edgeType)[edge] = (regions[i].layer == regions[j].layer) ? "direct-same-layer" : "direct-cross-layer";
                            directEdges++;
                        }
                        else if (result.connectionType == "blocked-clearable") {
                            (*edgeWeight)[edge] = result.blockedPixelsToRemove;
                            (*edgeCost)[edge] = result.pathCost;
                            (*edgeType)[edge] = "blocked-clearable";
                            blockedEdges++;
                        }
                        else if (result.connectionType == "cross-layer") {
                            (*edgeWeight)[edge] = 2;
                            (*edgeCost)[edge] = 2.0;
                            (*edgeType)[edge] = "cross-layer";
                            crossLayerEdges++;
                        }
                    }
                }
            }
        }
        
        cout << "边统计: 直接(" << directEdges << ") 可疏通(" << blockedEdges 
             << ") 跨层(" << crossLayerEdges << ") 总计(" << countEdges(graph) << ")" << endl;
    }
    
    // ===== 使用Lemon内置算法的全局分析功能 =====
    
    void analyzeLemonConnectivity() {
        cout << "\n=== Lemon内置算法连通性分析 ===" << endl;
        
        ListGraph::NodeMap<int> compMap(graph);
        int compCount = connectedComponents(graph, compMap);
        cout << "连通分量数: " << compCount << endl;
        
        cout << "图连通性: " << (connected(graph) ? "连通" : "非连通") << endl;
        
        if (countNodes(graph) > 1) {
            cout << "图是否双连通: " << (biNodeConnected(graph) ? "是" : "否") << endl;
            cout << "图是否双边连通: " << (biEdgeConnected(graph) ? "是" : "否") << endl;
        }
        
        analyzeComponentComposition(compMap, compCount);
        
        if (connected(graph)) {
            analyzeGraphMetrics();
        }
    }
    
    pair<vector<int>, double> findOptimalPath(int startRegionId, int endRegionId) {
        if (regionToNode.find(startRegionId) == regionToNode.end() ||
            regionToNode.find(endRegionId) == regionToNode.end()) {
            return {{}, -1.0};
        }
        
        ListGraph::Node startNode = regionToNode[startRegionId];
        ListGraph::Node endNode = regionToNode[endRegionId];
        
        Dijkstra<ListGraph, ListGraph::EdgeMap<double>> dijkstra(graph, *edgeCost);
        dijkstra.run(startNode);
        
        if (!dijkstra.reached(endNode)) {
            return {{}, -1.0};
        }
        
        vector<int> path;
        double totalCost = dijkstra.dist(endNode);
        
        ListGraph::Node current = endNode;
        while (current != INVALID) {
            path.push_back((*nodeRegionId)[current]);
            current = dijkstra.predNode(current);
        }
        
        reverse(path.begin(), path.end());
        return {path, totalCost};
    }
    
    vector<int> findMinHopPath(int startRegionId, int endRegionId) {
        if (regionToNode.find(startRegionId) == regionToNode.end() ||
            regionToNode.find(endRegionId) == regionToNode.end()) {
            return {};
        }
        
        ListGraph::Node startNode = regionToNode[startRegionId];
        ListGraph::Node endNode = regionToNode[endRegionId];
        
        Bfs<ListGraph> bfs(graph);
        bfs.run(startNode);
        
        if (!bfs.reached(endNode)) return {};
        
        vector<int> path;
        ListGraph::Node current = endNode;
        while (current != INVALID) {
            path.push_back((*nodeRegionId)[current]);
            current = bfs.predNode(current);
        }
        
        reverse(path.begin(), path.end());
        return path;
    }
    
    void findMinimumSpanningTree() {
        cout << "\n=== 最小生成树分析 ===" << endl;
        
        ListGraph::EdgeMap<bool> mstEdges(graph);
        double mstWeight = kruskal(graph, *edgeCost, mstEdges);
        
        int mstEdgeCount = 0;
        map<string, int> mstEdgeTypes;
        
        for (ListGraph::EdgeIt e(graph); e != INVALID; ++e) {
            if (mstEdges[e]) {
                mstEdgeCount++;
                mstEdgeTypes[(*edgeType)[e]]++;
            }
        }
        
        cout << "MST总权重: " << mstWeight << endl;
        cout << "MST边数: " << mstEdgeCount << endl;
        cout << "原图边数: " << countEdges(graph) << endl;
        
        if (countEdges(graph) > 0) {
            cout << "边减少率: " << (1.0 - (double)mstEdgeCount/countEdges(graph)) * 100 << "%" << endl;
        }
        
        cout << "MST边类型分布:" << endl;
        for (const auto& pair : mstEdgeTypes) {
            cout << "  " << pair.first << ": " << pair.second << endl;
        }
    }
    
    void performDFSAnalysis(int startRegionId) {
        if (regionToNode.find(startRegionId) == regionToNode.end()) {
            cout << "起始节点不存在" << endl;
            return;
        }
        
        cout << "\n=== DFS遍历分析 (从区域" << startRegionId << "开始) ===" << endl;
        
        ListGraph::Node startNode = regionToNode[startRegionId];
        
        Dfs<ListGraph> dfs(graph);
        dfs.run(startNode);
        
        vector<int> reachableNodes;
        for (ListGraph::NodeIt n(graph); n != INVALID; ++n) {
            if (dfs.reached(n)) {
                reachableNodes.push_back((*nodeRegionId)[n]);
            }
        }
        
        cout << "可达节点数: " << reachableNodes.size() << "/" << countNodes(graph) << endl;
        
        map<string, int> reachableByType;
        for (int regionId : reachableNodes) {
            ListGraph::Node node = regionToNode[regionId];
            reachableByType[(*nodeType)[node]]++;
        }
        
        cout << "可达节点类型分布:" << endl;
        for (const auto& pair : reachableByType) {
            cout << "  " << pair.first << ": " << pair.second << endl;
        }
    }
    
    // ===== 辅助分析函数 =====
    
    void analyzeComponentComposition(const ListGraph::NodeMap<int>& compMap, int compCount) {
        map<int, vector<int>> componentRegions;
        map<int, map<string, int>> componentTypeCount;
        
        for (ListGraph::NodeIt n(graph); n != INVALID; ++n) {
            int compId = compMap[n];
            int regionId = (*nodeRegionId)[n];
            string type = (*nodeType)[n];
            
            componentRegions[compId].push_back(regionId);
            componentTypeCount[compId][type]++;
        }
        
        cout << "\n连通分量详细分析:" << endl;
        for (int i = 0; i < compCount; i++) {
            cout << "分量 " << i << " (节点数: " << componentRegions[i].size() << "): ";
            
            bool first = true;
            for (const auto& pair : componentTypeCount[i]) {
                if (!first) cout << ", ";
                cout << pair.first << "(" << pair.second << ")";
                first = false;
            }
            cout << endl;
        }
    }
    
    void analyzeGraphMetrics() {
        cout << "\n图度量分析:" << endl;
        
        int totalDegree = 0;
        int maxDegree = 0;
        int minDegree = INT_MAX;
        
        for (ListGraph::NodeIt n(graph); n != INVALID; ++n) {
            int degree = 0;
            for (ListGraph::IncEdgeIt e(graph, n); e != INVALID; ++e) {
                degree++;
            }
            totalDegree += degree;
            maxDegree = max(maxDegree, degree);
            minDegree = min(minDegree, degree);
        }
        
        double avgDegree = (double)totalDegree / countNodes(graph);
        cout << "平均度数: " << avgDegree << endl;
        cout << "最大度数: " << maxDegree << endl;
        cout << "最小度数: " << minDegree << endl;
        
        int maxPossibleEdges = countNodes(graph) * (countNodes(graph) - 1) / 2;
        double density = (double)countEdges(graph) / maxPossibleEdges;
        cout << "图密度: " << density << endl;
    }
    
    // ===== 简化的连接分析函数 =====
    
    PathFindResult analyzeRegionConnection(const Region& regionA, const Region& regionB) {
        PathFindResult result;
        
        if (regionA.type != regionB.type) return result;
        
        int layerDiff = abs(regionA.layer - regionB.layer);
        
        if (layerDiff == 1) {
            int overlap = calculateCrossLayerOverlap(regionA, regionB);
            if (overlap >= overlapThreshold) {
                result.canConnect = true;
                result.connectionType = "cross-layer";
                result.pathCost = 2.0;
            }
            return result;
        }
        
        if (layerDiff > 1) return result;
        
        return analyzeSameLayerConnection(regionA, regionB);
    }
    
    PathFindResult analyzeSameLayerConnection(const Region& regionA, const Region& regionB) {
        PathFindResult result;
        
        int directOverlap = calculateSameLayerOverlap(regionA, regionB);
        if (directOverlap >= overlapThreshold) {
            result.canConnect = true;
            result.connectionType = "direct";
            result.pathCost = 1.0;
            return result;
        }
        
        // 简化的blocked路径检查
        auto centerA = regionA.getCenter();
        auto centerB = regionB.getCenter();
        double distance = sqrt(pow(centerA.first - centerB.first, 2) + pow(centerA.second - centerB.second, 2));
        
        if (distance < 15) { // 简单距离阈值
            result.canConnect = true;
            result.connectionType = "blocked-clearable";
            result.blockedPixelsToRemove = (int)(distance * 2);
            result.pathCost = distance * 0.5;
        }
        
        return result;
    }
    
    int calculateSameLayerOverlap(const Region& regionA, const Region& regionB) {
        set<pair<int, int>> pixelsA(regionA.pixels.begin(), regionA.pixels.end());
        int overlapCount = 0;
        for (const auto& pixel : regionB.pixels) {
            if (pixelsA.find(pixel) != pixelsA.end()) {
                overlapCount++;
            }
        }
        return overlapCount;
    }
    
    int calculateCrossLayerOverlap(const Region& regionA, const Region& regionB) {
        return calculateSameLayerOverlap(regionA, regionB);
    }
    
    void printGraphStats() {
        cout << "\n=== 图统计信息 ===" << endl;
        cout << "节点数: " << countNodes(graph) << endl;
        cout << "边数: " << countEdges(graph) << endl;
        
        map<string, int> edgeStats;
        double totalWeight = 0;
        
        for (ListGraph::EdgeIt e(graph); e != INVALID; ++e) {
            edgeStats[(*edgeType)[e]]++;
            totalWeight += (*edgeCost)[e];
        }
        
        cout << "平均边权重: " << (countEdges(graph) > 0 ? totalWeight / countEdges(graph) : 0) << endl;
        
        for (const auto& pair : edgeStats) {
            cout << pair.first << " 边数: " << pair.second << endl;
        }
    }
};

// ===== 测试用例创建 =====

void createTestPCBMatrix(vector<vector<vector<Cell>>>& matrix, int layers = 2, int width = 20, int height = 20) {
    matrix.resize(layers);
    
    for (int layer = 0; layer < layers; layer++) {
        matrix[layer].resize(width);
        for (int x = 0; x < width; x++) {
            matrix[layer][x].resize(height);
            for (int y = 0; y < height; y++) {
                matrix[layer][x][y] = Cell(layer, x, y);
            }
        }
    }
    
    // 创建测试区域
    auto gndRegion1 = make_shared<Region>(0, "GND", 0);
    auto gndRegion2 = make_shared<Region>(1, "GND", 0);
    auto gndRegion3 = make_shared<Region>(2, "GND", 1);
    auto gndRegion4 = make_shared<Region>(3, "GND", 0);
    
    // Layer 0 - 四个GND区域
    for (int x = 2; x <= 5; x++) {
        for (int y = 2; y <= 5; y++) {
            matrix[0][x][y].isEmpty = false;
            matrix[0][x][y].regionPtr = gndRegion1;
        }
    }
    
    for (int x = 12; x <= 15; x++) {
        for (int y = 12; y <= 15; y++) {
            matrix[0][x][y].isEmpty = false;
            matrix[0][x][y].regionPtr = gndRegion2;
        }
    }
    
    for (int x = 2; x <= 5; x++) {
        for (int y = 12; y <= 15; y++) {
            matrix[0][x][y].isEmpty = false;
            matrix[0][x][y].regionPtr = gndRegion4;
        }
    }
    
    // Layer 1 - 与Layer 0有重叠的GND区域
    for (int x = 3; x <= 6; x++) {
        for (int y = 3; y <= 6; y++) {
            matrix[1][x][y].isEmpty = false;
            matrix[1][x][y].regionPtr = gndRegion3;
        }
    }
    
    // 添加blocked区域
    for (int x = 6; x <= 11; x++) {
        for (int y = 7; y <= 9; y++) {
            matrix[0][x][y].isEmpty = false;
            matrix[0][x][y].isBlocked = true;
        }
    }
    
    cout << "创建了测试PCB矩阵 (" << layers << "层, " << width << "x" << height << ")" << endl;
}

// ===== 主函数 =====

int main() {
    cout << "=== PCB连通性图构建 - 修复版 ===" << endl;
    
    // 创建测试数据
    vector<vector<vector<Cell>>> testMatrix;
    createTestPCBMatrix(testMatrix, 2, 20, 20);
    
    // 初始化图构建器
    PCBConnectivityGraphBuilder builder(3, 30);
    
    // 构建图
    builder.extractRegionsFromMatrix(testMatrix);
    builder.buildConnectivityGraph();
    
    // 使用Lemon内置算法进行全面分析
    builder.analyzeLemonConnectivity();
    builder.findMinimumSpanningTree();
    builder.performDFSAnalysis(0);
    
    // 路径查找测试
    cout << "\n=== 路径查找测试 ===" << endl;
    
    auto optimalPath = builder.findOptimalPath(0, 1);
    if (!optimalPath.first.empty()) {
        cout << "最优路径 (区域0->1): ";
        for (int i = 0; i < optimalPath.first.size(); i++) {
            if (i > 0) cout << " -> ";
            cout << optimalPath.first[i];
        }
        cout << " (成本: " << optimalPath.second << ")" << endl;
    }
    
    auto minHopPath = builder.findMinHopPath(0, 1);
    if (!minHopPath.empty()) {
        cout << "最少跳数路径 (区域0->1): ";
        for (int i = 0; i < minHopPath.size(); i++) {
            if (i > 0) cout << " -> ";
            cout << minHopPath[i];
        }
        cout << " (跳数: " << minHopPath.size() - 1 << ")" << endl;
    }
    
    // 打印图统计信息
    builder.printGraphStats();
    
    return 0;
}