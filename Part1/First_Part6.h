#pragma once
#include <string>
#include <vector>
#include <memory>
#include <iosfwd>
#include <cstddef>

// ===== AST 节点 =====
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

// ===== 板层信息 (layers ...) =====
struct LayerInfo {

    int         id = -1;          //  0, 31, 44 ...
    std::string name;             //  F.Cu, B.Cu, Edge.Cuts ...
    std::string kind;             //  signal, user
    std::string description;      //  "Top", "Bottom", "B.Courtyard"Ϊ
};

class KiCadParser {

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
    //segment
    void findAndNumberSegments();
    void findAndNumberVias();
    void findAndNumberFootprints();

public:
    KiCadParser();

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


};
#pragma once
