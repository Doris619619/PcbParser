#ifndef FIRST_PART12_H
#define FIRST_PART12_H

#include <string>
#include <vector>
#include <memory>
#include <iosfwd>
#include <cstddef>

// ===== AST �ڵ� =====
struct Node {

    std::string name;
    std::vector<std::string> parameters;
    std::vector<std::shared_ptr<Node>> children;

    Node(const std::string& n = "") : name(n) {}

    void addParameter(const std::string& param) {
        parameters.push_back(param);
    }

    void addChild(std::shared_ptr<Node> child) {
        children.push_back(child);
    }
};

// ===== �����Ϣ (layers ...) =====
struct LayerInfo {

    int         id = -1;          //  0, 31, 44 ...
    std::string name;             //  F.Cu, B.Cu, Edge.Cuts ...
    std::string kind;             //  signal, user
    std::string description;      //  "Top", "Bottom", "B.Courtyard"?
};

struct NetInfo {
    int         id = -1;      // e.g. 0,1,2...
    std::string name;         // e.g. "VIN","GND"�������ǿմ�
};

class KiCadParser {

public:
    KiCadParser();

    // [NEW] �����Ӿ��ν��
    struct BoardBBox {
        bool   valid = false;
        double minX = 0.0, minY = 0.0, maxX = 0.0, maxY = 0.0;
        double width()  const { return valid ? (maxX - minX) : 0.0; }
        double height() const { return valid ? (maxY - minY) : 0.0; }
    };

    // �򵥶�ά��
    struct Point2D { double x = 0.0; double y = 0.0; };

    // ȡ���ΰ����ĸ��սǣ��������¡����ϡ����ϡ����£�
    bool getBoardCorners(Point2D& bottomLeft,
        Point2D& topLeft,
        Point2D& topRight,
        Point2D& bottomRight) const;

    // === Net �ӿڣ�������===
    const std::vector<NetInfo>& getBoardNets() const;
    void printBoardNets() const;

    // [NEW] ���� Edge.Cuts ������Ӿ���
    BoardBBox getBoardBBox() const;

    // [NEW] ֱ���ÿ�ߣ���ʧ�ܷ��� false��
    bool getBoardSize(double& width, double& height) const;

    bool saveAsKicadPcb(const std::string& outPath,
        bool stripTempIds = true,
        int indentStep = 2);
    size_t removeViaByTstamp(const std::string& t);

    size_t removeSegmentByTstamp(const std::string& t);

    const std::vector<LayerInfo>& getBoardLayers() const;
    void printBoardLayers() const;
    int nextSegmentId() const;

    int addSegment(double startX, double startY, double endX, double endY,
        double width, const std::string& layer, int net, const std::string& tstamp);

    int addSegmentSimple(double startX, double startY, double endX, double endY,
        double width, const std::string& layer,
        int net, const std::string& tstamp);

    int addVia(double x, double y,
        double size, double drill,
        const std::vector<std::string>& layers,
        bool isFree,
        int net,
        const std::string& tstamp);

    int addViaSimple(double x, double y,
        double size, double drill,
        const std::string& layer1,
        const std::string& layer2_or_dash,
        int freeFlag,
        int net,
        const std::string& tstamp);

    void stripAllTempIds();
    bool parseFile(const std::string& filename);
    const std::vector<std::shared_ptr<Node>>& getAllSegments() const;
    const std::vector<std::shared_ptr<Node>>& getAllVias() const;
    const std::vector<std::shared_ptr<Node>>& getAllFootprints() const;
    const std::vector<std::shared_ptr<Node>>& getAllNodes() const;
    size_t removeViaById(int id);
    size_t removeSegmentById(int id);
    void printSegmentInfo(const std::shared_ptr<Node>& segment);
    void printViaInfo(const std::shared_ptr<Node>& via);
    void printFootprintInfo(const std::shared_ptr<Node>& footprint);
    void printAllSegments();
    void printAllVias();
    void printAllFootprints();
    void printStructure(const std::shared_ptr<Node>& node = nullptr, int depth = 0, int maxDepth = 3);

private:
    std::shared_ptr<Node> root;
    size_t pos;
    std::string content;
    int segmentCounter;
    int viaCounter;
    std::vector<std::shared_ptr<Node>> allSegments;
    std::vector<std::shared_ptr<Node>> allVias;
    std::vector<std::shared_ptr<Node>> allFootprints;
    std::vector<std::shared_ptr<Node>> allNodes;
    std::vector<LayerInfo> boardLayers;
    std::vector<NetInfo> boardNets;

    void parseBoardNets();
    static std::string unquote(const std::string& s);
    static bool isEdgeCuts(const std::shared_ptr<Node>& n);
    static bool readXY(const std::shared_ptr<Node>& n, const char* childName, double& x, double& y);
    static inline void pushPoint(double x, double y, BoardBBox& bbox);
    size_t removeNodesByTstampRec(const std::shared_ptr<Node>& cur,
        const std::string& tstamp,
        const char* expectName);
    void parseBoardLayers();
    static std::string fmtNum(double v);
    static int extractTagId(const std::string& p, const char* kind);
    int nextViaId() const;
    static inline bool hasExactTag(const std::shared_ptr<Node>& n, const char* kind, int id);
    size_t removeNodesByTagRec(const std::shared_ptr<Node>& cur,
        const char* kind,
        int id,
        const char* expectName);
    void rebuildCaches();
    static bool isNumberToken(const std::string& s);
    static std::string escapeForQuotes(const std::string& s);
    static bool shouldQuote(const std::string& s);
    static std::string renderParam(const std::string& p);
    void writeNodeRec(std::ostream& os,
        const std::shared_ptr<Node>& n,
        int indent,
        int indentStep);
    void stripAllTempIdsRec(const std::shared_ptr<Node>& n);
    void skipWhitespace();
    std::string readBareToken();
    char currentChar();
    char nextChar();
    std::string readQuotedString();
    std::string readToken();
    std::shared_ptr<Node> parseNode();
    void findNodesRecursive(const std::shared_ptr<Node>& node,
        const std::string& targetName,
        std::vector<std::shared_ptr<Node>>& results);
    void collectAllNodes(const std::shared_ptr<Node>& node);
    void findAndNumberSegments();
    void findAndNumberVias();
    void findAndNumberFootprints();
};

#endif // FIRST_PART12_H#pragma once
