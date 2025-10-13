#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include <cctype>
#include <algorithm>
#include <functional>
#include <unordered_set>

#include <iomanip>
#include <climits>
#include "First_Part12.h"


//�ڱ�
// === Tag helpers: [[FPID:x]] / [[SEID:x]] / [[VIID:x]] ===
static inline std::string makeTag(const char* kind, int id) {
    // kind = "FPID" / "SEID" / "VIID"
    return std::string("[[") + kind + ":" + std::to_string(id) + "]]";
}
static inline bool isTagged(const std::string& s, const char* kind) {
    const std::string p = std::string("[[") + kind + ":";
    return s.rfind(p, 0) == 0 && s.size() >= p.size() + 2 && s.substr(s.size() - 2) == "]]";
}


// ===== ɾ���������� tstamp ��ȷƥ�� =====
// ��ɾ�� name==expectName �����ӽڵ���� (tstamp <exact>) �Ľڵ�
size_t KiCadParser::removeNodesByTstampRec(const std::shared_ptr<Node>& cur,
    const std::string& tstamp,
    const char* expectName) {
    if (!cur) return 0;
    size_t removed = 0;
    auto& ch = cur->children;
    const size_t before = ch.size();
    ch.erase(std::remove_if(ch.begin(), ch.end(),
        [&](const std::shared_ptr<Node>& c) {
            if (!c) return false;
            if (expectName && c->name != expectName) return false;
            for (const auto& gc : c->children) {
                if (gc && gc->name == "tstamp" && !gc->parameters.empty() && gc->parameters[0] == tstamp)
                    return true;
            }
            return false;
        }), ch.end());
    removed += (before - ch.size());
    for (auto& c : ch) removed += removeNodesByTstampRec(c, tstamp, expectName);
    return removed;
}

// === Net: ����������� (net <id> "<name>") ===
void KiCadParser::parseBoardNets() {
    boardNets.clear();
    if (!root) return;

    for (const auto& ch : root->children) {
        if (!ch || ch->name != "net") continue;

        NetInfo info;
        // (net 21 "GND")
        if (!ch->parameters.empty()) {
            try { info.id = std::stoi(ch->parameters[0]); }
            catch (...) { info.id = -1; }
        }
        if (ch->parameters.size() >= 2) {
            info.name = ch->parameters[1]; // �����ǿմ� ""
        }
        boardNets.push_back(std::move(info));
    }
}

const std::vector<NetInfo>& KiCadParser::getBoardNets() const {
    return boardNets;
}

void KiCadParser::printBoardNets() const {
    for (const auto& n : boardNets) {
        std::cout << "Net " << n.id << " : \"" << n.name << "\"\n";
    }
}

// �������� (layers ...)��ֻɨ root ��ֱ���ӽڵ㣬������ pad �ڵ� (layers ...) ����
void KiCadParser::parseBoardLayers() {
    boardLayers.clear();
    if (!root) return;

    // ֻ�� root->children ��������Ϊ "layers" �����ӽڵ�Ŀ飨����������������̬��
    for (const auto& ch : root->children) {
        if (!ch) continue;
        if (ch->name != "layers") continue;
        if (ch->children.empty()) continue; // pad �� (layers "F.Cu" ...) û���ӽڵ�

        for (const auto& lay : ch->children) {
            if (!lay) continue;

            // �������磺 (0 "F.Cu" signal "Top")
            // ��������name=��ţ�parameters[0]=������parameters[1]=���ͣ�parameters[2]=��ѡ����
            int id = -1;
            try { id = std::stoi(lay->name); }
            catch (...) { /* �����־����� */ continue; }

            LayerInfo info;
            info.id = id;
            if (!lay->parameters.empty())      info.name = lay->parameters[0];
            if (lay->parameters.size() >= 2)   info.kind = lay->parameters[1];
            if (lay->parameters.size() >= 3)   info.description = lay->parameters[2];

            boardLayers.push_back(std::move(info));
        }
        break; // ���� (layers ...) ֻ����һ�����ҵ���ͣ
    }
}



//���via����
// ���� ����ת�ַ�����������ԭ�����ܶ��������ӣ�����ѧ�����������ʾ��ȣ�
std::string KiCadParser::fmtNum(double v) {
    std::ostringstream oss;
    oss << std::setprecision(15) << std::noshowpoint << std::defaultfloat << v;
    return oss.str();
}

// ���� �� [[KIND:123]] ����ȡ id����ȡʧ�ܷ��� -1
int KiCadParser::extractTagId(const std::string& p, const char* kind) {
    if (!isTagged(p, kind)) return -1;
    const size_t colon = p.find(':');
    const size_t end = p.rfind("]]");
    if (colon == std::string::npos || end == std::string::npos || end <= colon + 1) return -1;
    try { return std::stoi(p.substr(colon + 1, end - colon - 1)); }
    catch (...) { return -1; }
}

// ���� ������һö via �ı�ţ����� [[VIID:x]] �����ֵ + 1����û�У���� 1 ��ʼ
int KiCadParser::nextViaId() const {
    int mx = 0;
    for (const auto& v : allVias) {
        if (!v || v->parameters.empty()) continue;
        int id = extractTagId(v->parameters.front(), "VIID");
        if (id > mx) mx = id;
    }
    return mx + 1;
}


// === ɾ����أ�ƥ����ݹ��Ƴ� ===

// �жϽڵ���ײ����Ƿ��������ľ�ȷ�ڱ������� [[VIID:12]]
inline bool KiCadParser::hasExactTag(const std::shared_ptr<Node>& n, const char* kind, int id) {
    if (!n || n->parameters.empty()) return false;
    return n->parameters.front() == makeTag(kind, id);
}

// �ݹ����������ɾ������������ (expectName && �ײ��� == [[KIND:id]]) �Ľڵ�
// ����ɾ���������������� 0 / 1 / �����
size_t KiCadParser::removeNodesByTagRec(const std::shared_ptr<Node>& cur,
    const char* kind,
    int id,
    const char* expectName) {
    if (!cur) return 0;
    size_t removed = 0;

    // ���ڵ�ǰ�ڵ��ֱ�� children ��ɾ
    auto& ch = cur->children;
    const size_t before = ch.size();
    ch.erase(std::remove_if(ch.begin(), ch.end(),
        [&](const std::shared_ptr<Node>& c) {
            if (!c) return false;
            if (expectName && c->name != expectName) return false;
            return hasExactTag(c, kind, id);
        }),
        ch.end());
    removed += (before - ch.size());

    // �ٵݹ鴦��ʣ�ຢ��
    for (auto& c : ch) {
        removed += removeNodesByTagRec(c, kind, id, expectName);
    }
    return removed;
}

// ɾ�������±�ţ�ֻ���ؽ����桱�б�����������Ų���
void KiCadParser::rebuildCaches() {
    allNodes.clear();
    collectAllNodes(root);

    allVias.clear();
    findNodesRecursive(root, "via", allVias);

    allSegments.clear();
    findNodesRecursive(root, "segment", allSegments);

    allFootprints.clear();
    findNodesRecursive(root, "footprint", allFootprints);
}


// ----- ���л��������ж����� -----
bool KiCadParser::isNumberToken(const std::string& s) {
    if (s.empty()) return false;
    bool seenDigit = false, seenDot = false, seenExp = false;
    size_t i = 0; if (s[i] == '+' || s[i] == '-') ++i;
    for (; i < s.size(); ++i) {
        char c = s[i];
        if (std::isdigit((unsigned char)c)) { seenDigit = true; continue; }
        if (c == '.') { if (seenDot || seenExp) return false; seenDot = true; continue; }
        if (c == 'e' || c == 'E') {
            if (seenExp || !seenDigit) return false;
            seenExp = true; seenDigit = false;
            if (i + 1 < s.size() && (s[i + 1] == '+' || s[i + 1] == '-')) ++i;
            continue;
        }
        return false;
    }
    return seenDigit;
}

// ----- ת�嵽˫���� -----
std::string KiCadParser::escapeForQuotes(const std::string& s) {
    std::string out; out.reserve(s.size() + 8);
    for (char c : s) { if (c == '"' || c == '\\') out.push_back('\\'); out.push_back(c); }
    return out;
}

// ----- �Ƿ���Ҫ���� -----
// ���򣺿մ� �� �������ţ����� �� �����ţ�
// ֻ�� [a-z0-9_\-./+:] �� �����ţ�ȫ��Сд/��ȫ�ַ�����
// ����д��հ�/�������� �� ���š�
bool KiCadParser::shouldQuote(const std::string& s) {
    if (s.empty()) return true;           // �ؼ������ ""������ group "" / net_name "" / outputdirectory ""
    if (isNumberToken(s)) return false;
    bool hasUpper = false;
    for (unsigned char uc : s) {
        char c = (char)uc;
        if (std::isupper(uc)) hasUpper = true;
        if (!(std::islower(uc) || std::isdigit(uc) ||
            c == '_' || c == '-' || c == '.' || c == '/' || c == '+' || c == ':')) {
            return true; // �пհס����š���б�ܻ���������
        }
    }
    return hasUpper; // �д�дҲ�����ţ��� "B.Cu", "Top"��
}

// ----- ������Ⱦ -----
std::string KiCadParser::renderParam(const std::string& p) {
    if (p.empty()) return "\"\"";                      // �� ���ַ�����ǿ����� ""
    if (!shouldQuote(p)) return p;
    return std::string("\"") + escapeForQuotes(p) + "\"";
}

// ----- �ݹ�д�ڵ㣨���� hide ��λ�ã� -----
void KiCadParser::writeNodeRec(std::ostream& os,
    const std::shared_ptr<Node>& n,
    int indent,
    int indentStep) {
    if (!n) return;

    std::string indentStr((size_t)indent, ' ');

    // 1) �ȰѲ������ "hide" �����
    std::vector<std::string> paramsNoHide;
    int hideCount = 0;
    for (const auto& p : n->parameters) {
        if (p == "hide") ++hideCount;
        else paramsNoHide.push_back(p);
    }

    // 2) д��(name + ��ͨ����)��
    os << indentStr << "(" << n->name;
    for (const auto& p : paramsNoHide)
        os << " " << renderParam(p);

    // 3) ���û���ӽڵ㣺�� hide �Է�����β�����ּ���
    if (n->children.empty()) {
        while (hideCount-- > 0) os << " hide";
        os << ")\n";
        return;
    }

    // 4) ���ӽڵ㣺���У�Ȼ����д���һ���ӽڵ�����һ�� "hide"
    os << "\n";

    bool injectedHide = false;
    for (size_t i = 0; i < n->children.size(); ++i) {
        const auto& ch = n->children[i];
        writeNodeRec(os, ch, indent + indentStep, indentStep);

        if (!injectedHide && hideCount > 0) {
            std::string ind2((size_t)(indent + indentStep), ' ');
            while (hideCount-- > 0) {
                os << ind2 << "hide\n";  // �����ڵ�һ�����������֮��
            }
            injectedHide = true;
        }
    }

    os << indentStr << ")\n";
}







//�ڱ�
// --- �� KiCadParser::private: �� ---
void KiCadParser::stripAllTempIdsRec(const std::shared_ptr<Node>& n) {
    if (!n) return;
    if (!n->parameters.empty()) {
        const std::string& p0 = n->parameters.front();
        if (isTagged(p0, "FPID") || isTagged(p0, "SEID") || isTagged(p0, "VIID")) {
            const_cast<std::vector<std::string>&>(n->parameters).erase(n->parameters.begin());
        }
    }
    for (auto& ch : n->children) stripAllTempIdsRec(ch);
}




/*
void skipWhitespace() {
    while (pos < content.size() && std::isspace(content[pos])) {
        pos++;
    }
}
*/

// 1) ����� isspace �Ĳ�����ת�� unsigned char����ֹ����/��λ�ֽ�����UB
void KiCadParser::skipWhitespace() {
    while (pos < content.size() && std::isspace(static_cast<unsigned char>(content[pos]))) {
        ++pos;
    }
}

// 2) �������ӵ�ǰλ�ö�ȡһ������ token������������ͷ�հף���������
std::string KiCadParser::readBareToken() {
    std::string token;
    while (pos < content.size()) {
        char c = currentChar();
        if (std::isspace(static_cast<unsigned char>(c)) || c == '(' || c == ')') break;
        token += c;
        ++pos;
    }
    return token;
}

char KiCadParser::currentChar() {
    if (pos < content.size()) {
        return content[pos];
    }
    return '\0';
}

char KiCadParser::nextChar() {
    pos++;
    return currentChar();
}

std::string KiCadParser::readQuotedString() {
    if (currentChar() != '"') return "";
    ++pos;                       // �Ե���ͷ�� "
    std::string out;
    while (pos < content.size() && currentChar() != '"') {
        if (currentChar() == '\\') ++pos;     // ����ת��
        if (pos < content.size()) {
            out.push_back(currentChar());
            ++pos;
        }
    }
    if (currentChar() == '"') ++pos;          // �Ե���β�� "
    return out;                               // �����ǿմ� ""
}



std::string KiCadParser::readToken() {
    skipWhitespace();

    if (pos >= content.size()) {
        return "";
    }

    char c = currentChar();
    if (c == '(' || c == ')') {
        pos++;
        return std::string(1, c);
    }

    if (c == '"') {
        return "\"" + readQuotedString() + "\"";
    }

    std::string token;
    while (pos < content.size()) {
        c = currentChar();
        if (std::isspace(c) || c == '(' || c == ')') {
            break;
        }
        token += c;
        pos++;
    }

    return token;
}

std::shared_ptr<Node> KiCadParser::parseNode() {
    skipWhitespace();

    if (currentChar() != '(') {
        return nullptr;
    }

    ++pos;
    skipWhitespace();

    std::string name;
    while (pos < content.size()) {
        char c = currentChar();
        if (std::isspace(static_cast<unsigned char>(c)) || c == '(' || c == ')') break;
        name += c;
        ++pos;
    }
    if (name.empty()) return nullptr;

    auto node = std::make_shared<Node>(name);

    while (pos < content.size()) {
        skipWhitespace();

        if (currentChar() == ')') { ++pos; break; }

        if (currentChar() == '(') {
            auto child = parseNode();
            if (child) node->addChild(child);
        }
        else {
            if (currentChar() == '"') {
                // �� �����մ�����
                std::string q = readQuotedString();
                node->addParameter(q);
            }
            else {
                std::string tok = readBareToken();   // ��û������������� readToken()
                if (!tok.empty()) node->addParameter(tok);
            }
        }
    }

    return node;
}


/*
void findNodesRecursive(const std::shared_ptr<Node>& node,
    const std::string& targetName,
    std::vector<std::shared_ptr<Node>>& results) {
    if (!node) return;

    if (node->name == targetName) {
        results.push_back(node);
    }

    for (const auto& child : node->children) {
        findNodesRecursive(child, targetName, results);
    }
}
*/
/*
void findNodesRecursive(const std::shared_ptr<Node>& node,
    const std::string& targetName,
    std::vector<std::shared_ptr<Node>>& results) {
    if (!node) return;

    // �����нڵ���ӵ� allNodes ��
    allNodes.push_back(node);

    // �����ض����ƵĽڵ㲢������ӵ�����б���
    if (node->name == targetName) {
        results.push_back(node);
    }

    // �ݹ�����ӽڵ�
    for (const auto& child : node->children) {
        findNodesRecursive(child, targetName, results);
    }
}
*/

// ֻ�����������ҡ��������� allNodes �����κζ���
void KiCadParser::findNodesRecursive(const std::shared_ptr<Node>& node,
    const std::string& targetName,
    std::vector<std::shared_ptr<Node>>& results) {
    if (!node) return;
    if (node->name == targetName) {
        results.push_back(node);
    }
    for (const auto& child : node->children) {
        findNodesRecursive(child, targetName, results);
    }
}

// ר�������ռ����нڵ�
void KiCadParser::collectAllNodes(const std::shared_ptr<Node>& node) {
    if (!node) return;
    allNodes.push_back(node);
    for (const auto& child : node->children) {
        collectAllNodes(child);
    }
}




//������segment���
void KiCadParser::findAndNumberSegments() {
    segmentCounter = 0;
    allSegments.clear();
    findNodesRecursive(root, "segment", allSegments);

    for (auto& segment : allSegments) {
        segmentCounter++;
        std::string tag = makeTag("SEID", segmentCounter);
        if (segment->parameters.empty() || !isTagged(segment->parameters.front(), "SEID")) {
            segment->parameters.insert(segment->parameters.begin(), tag);
        }
        else {
            segment->parameters.front() = tag; // �����ڱ��򸲸�
        }

    }
}



void KiCadParser::findAndNumberVias() {
    viaCounter = 0;
    allVias.clear();
    findNodesRecursive(root, "via", allVias);
    for (auto& via : allVias) {
        ++viaCounter;
        std::string tag = makeTag("VIID", viaCounter);
        if (via->parameters.empty() || !isTagged(via->parameters.front(), "VIID")) {
            via->parameters.insert(via->parameters.begin(), tag);
        }
        else {
            via->parameters.front() = tag;
        }

    }
}


void KiCadParser::findAndNumberFootprints() {
    int footprintCounter = 0;  // ������footprint���
    allFootprints.clear();  // ���ԭ�е�footprints
    findNodesRecursive(root, "footprint", allFootprints);  // ��������footprint�ڵ�

    for (auto& footprint : allFootprints) {
        footprintCounter++;
        std::string tag = makeTag("FPID", footprintCounter);
        if (footprint->parameters.empty() || !isTagged(footprint->parameters.front(), "FPID")) {
            footprint->parameters.insert(footprint->parameters.begin(), tag);
        }
        else {
            footprint->parameters.front() = tag;
        }

    }
}

KiCadParser::KiCadParser() : pos(0), segmentCounter(0), viaCounter(0) {}

// ===== ���� .kicad_pcb =====
// Ĭ�ϱ���ǰ��� [[FPID:]] / [[SEID:]] / [[VIID:]] �ڱ���������Ⱦ�ļ�
bool KiCadParser::saveAsKicadPcb(const std::string& outPath,
    bool stripTempIds,
    int indentStep) {
    if (!root) return false;
    if (stripTempIds) stripAllTempIds();
    std::ofstream ofs(outPath, std::ios::binary);
    if (!ofs.is_open()) {
        std::cerr << "�޷�д���ļ�: " << outPath << "\n";
        return false;
    }
    writeNodeRec(ofs, root, /*indent=*/0, indentStep);
    return true;
}

// ===== �� tstamp ɾ��������δ��Ŷ��� =====
size_t KiCadParser::removeViaByTstamp(const std::string& t) {
    size_t n = removeNodesByTstampRec(root, t, "via");
    if (n > 0) rebuildCaches();
    return n;
}
size_t KiCadParser::removeSegmentByTstamp(const std::string& t) {
    size_t n = removeNodesByTstampRec(root, t, "segment");
    if (n > 0) rebuildCaches();
    return n;
}

//�洢layer vector
// ��ȡ�������б�
const std::vector<LayerInfo>& KiCadParser::getBoardLayers() const {
    return boardLayers;
}

// ���Դ�ӡ�����һ��
void KiCadParser::printBoardLayers() const {
    std::cout << "\n=== ��㶨�� (�� " << boardLayers.size() << " ��) ===\n";
    for (const auto& L : boardLayers) {
        std::cout << std::setw(2) << L.id << "  "
            << L.name << "  "
            << L.kind;
        if (!L.description.empty()) {
            std::cout << "  \"" << L.description << "\"";
        }
        std::cout << "\n";
    }
}


//����segment

// ������һö segment �ı�ţ����� [[SEID:x]] �����ֵ + 1����û�У���� 1 ��ʼ
int KiCadParser::nextSegmentId() const {
    int mx = 0;
    for (const auto& s : allSegments) {
        if (!s || s->parameters.empty()) continue;
        int id = extractTagId(s->parameters.front(), "SEID");
        if (id > mx) mx = id;
    }
    return mx + 1;
}

// ͨ�ð棺layers �ɴ� 1 ����ͨ�� 1 ������������ segment �ı�ţ�SEID��
int KiCadParser::addSegment(double startX, double startY, double endX, double endY,
    double width, const std::string& layer, int net, const std::string& tstamp)
{
    if (!root) return -1;

    // 1) �����±�����ڱ�
    const int newId = nextSegmentId();
    const std::string tag = makeTag("SEID", newId);

    // 2) ���� segment �ڵ㣨�ṹ��ԭ�ļ�һ�£�segment ������ڱ����������඼���ӽڵ㣩
    auto segment = std::make_shared<Node>("segment");
    segment->parameters.push_back(tag);   // [[SEID:x]]

    // (start x y)
    {
        auto startN = std::make_shared<Node>("start");
        startN->parameters.push_back(fmtNum(startX));
        startN->parameters.push_back(fmtNum(startY));
        segment->addChild(startN);
    }

    // (end x y)
    {
        auto endN = std::make_shared<Node>("end");
        endN->parameters.push_back(fmtNum(endX));
        endN->parameters.push_back(fmtNum(endY));
        segment->addChild(endN);
    }

    // (width w)
    {
        auto widthN = std::make_shared<Node>("width");
        widthN->parameters.push_back(fmtNum(width));
        segment->addChild(widthN);
    }

    // (layer)
    {
        auto layerN = std::make_shared<Node>("layer");
        layerN->parameters.push_back(layer);
        segment->addChild(layerN);
    }

    // (net n)
    {
        auto netN = std::make_shared<Node>("net");
        netN->parameters.push_back(std::to_string(net));
        segment->addChild(netN);
    }

    // (tstamp xxx)
    {
        auto tN = std::make_shared<Node>("tstamp");
        tN->parameters.push_back(tstamp);
        segment->addChild(tN);
    }

    // 3) �嵽���㣺��ԭ�� segment һ������ root ֱ��
    root->addChild(segment);

    // 4) ˢ�»��棨���ر�ţ�������������� ID �ȶ���
    rebuildCaches();

    return newId;
}

// ��ݰ棺������ʽ���루��(layer)���մ���-����ʾһ���㣩
int KiCadParser::addSegmentSimple(double startX, double startY, double endX, double endY,
    double width, const std::string& layer,
    int net, const std::string& tstamp)
{
    return addSegment(startX, startY, endX, endY, width, layer, net, tstamp);
}



// ͨ�ð棺layers �ɴ� 0/1/N ����ͨ�� 1 �� 2 ������isFree= true/false
// ������ via �ı�ţ�VIID��
int KiCadParser::addVia(double x, double y,
    double size, double drill,
    const std::vector<std::string>& layers,
    bool isFree,
    int net,
    const std::string& tstamp)
{
    if (!root) return -1;

    // 1) �����±�����ڱ�
    const int newId = nextViaId();
    const std::string tag = makeTag("VIID", newId);

    // 2) ���� via �ڵ㣨�ṹ��ԭ�ļ�һ�£�via ������ڱ����������඼���ӽڵ㣩
    auto via = std::make_shared<Node>("via");
    via->parameters.push_back(tag);   // [[VIID:x]]

    // (at x y)
    {
        auto atN = std::make_shared<Node>("at");
        atN->parameters.push_back(fmtNum(x));
        atN->parameters.push_back(fmtNum(y));
        via->addChild(atN);
    }
    // (size s)
    {
        auto sN = std::make_shared<Node>("size");
        sN->parameters.push_back(fmtNum(size));
        via->addChild(sN);
    }
    // (drill d)
    {
        auto dN = std::make_shared<Node>("drill");
        dN->parameters.push_back(fmtNum(drill));
        via->addChild(dN);
    }
    // (layers ...) ���� ��Щֻ�� 1 ����������Щ 2 �������ǰ����� vector ��̬����
    if (!layers.empty()) {
        auto layN = std::make_shared<Node>("layers");
        for (const auto& L : layers) {
            if (!L.empty()) layN->parameters.push_back(L); // �ڲ��治������
        }
        // ����ȷʵ�в���ʱ�Ų���
        if (!layN->parameters.empty()) via->addChild(layN);
    }
    // (free) ��ѡ
    if (isFree) {
        via->addChild(std::make_shared<Node>("free"));  // �޲���
    }
    // (net n)
    {
        auto nN = std::make_shared<Node>("net");
        nN->parameters.push_back(std::to_string(net));
        via->addChild(nN);
    }
    // (tstamp xxx)
    {
        auto tN = std::make_shared<Node>("tstamp");
        tN->parameters.push_back(tstamp);
        via->addChild(tN);
    }

    // 3) �嵽���㣺��ԭ�� via һ������ root ֱ��
    root->addChild(via);

    // 4) ˢ�»��棨���ر�ţ�������������� ID �ȶ���
    rebuildCaches();

    return newId;
}

// ��ݰ棺����˵�ġ�������ʽ��������layer2 ���մ��� "-" ����ֻ�� 1 �� layer��freeFlag 0/1��
int KiCadParser::addViaSimple(double x, double y,
    double size, double drill,
    const std::string& layer1,
    const std::string& layer2_or_dash,
    int freeFlag,
    int net,
    const std::string& tstamp)
{
    std::vector<std::string> layers;
    if (!layer1.empty()) layers.push_back(layer1);
    if (!layer2_or_dash.empty() && layer2_or_dash != "-") layers.push_back(layer2_or_dash);
    return addVia(x, y, size, drill, layers, freeFlag != 0, net, tstamp);
}

void KiCadParser::stripAllTempIds() { stripAllTempIdsRec(root); }

bool KiCadParser::parseFile(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "�޷����ļ�: " << filename << std::endl;
        return false;
    }

    content.assign((std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>());
    file.close();

    pos = 0;
    root = parseNode();

    if (root) {
        allNodes.clear();                 // �����
        collectAllNodes(root);            // һ�����ռ�ȫ���ڵ�

        parseBoardLayers();
        parseBoardNets();


        findAndNumberSegments();
        findAndNumberVias();
        findAndNumberFootprints();
    }

    return root != nullptr;
}

// ��ȡ����segment
const std::vector<std::shared_ptr<Node>>& KiCadParser::getAllSegments() const {
    return allSegments;
}

// ��ȡ����via
const std::vector<std::shared_ptr<Node>>& KiCadParser::getAllVias() const {
    return allVias;
}

// ��ȡ����footprint
const std::vector<std::shared_ptr<Node>>& KiCadParser::getAllFootprints() const {
    return allFootprints;
}

//��ȡ���нڵ�Node
const std::vector<std::shared_ptr<Node>>& KiCadParser::getAllNodes() const {
    return allNodes;
}

// === ����ɾ�� API ===
// ���ݱ��ɾ������ƥ��� via������ [[VIID:x]]��������ɾ������
size_t KiCadParser::removeViaById(int id) {
    // ע�⣺ɾ�������ڱ����ڣ��벻Ҫ�ڵ���ǰ�� stripAllTempIds()
    size_t n = removeNodesByTagRec(root, "VIID", id, "via");
    if (n > 0) rebuildCaches();
    return n;
}

// ���ݱ��ɾ������ƥ��� segment������ [[SEID:x]]��������ɾ������
size_t KiCadParser::removeSegmentById(int id) {
    size_t n = removeNodesByTagRec(root, "SEID", id, "segment");
    if (n > 0) rebuildCaches();
    return n;
}


// ��ӡsegment��Ϣ
void KiCadParser::printSegmentInfo(const std::shared_ptr<Node>& segment) {
    if (!segment || segment->name != "segment") {
        std::cout << "������Ч��segment�ڵ�" << std::endl;
        return;
    }

    std::cout << "Segment #" << segment->parameters[0] << ":" << std::endl;

    if (segment->parameters.size() > 1) {
        std::cout << "  ����: ";
        for (size_t i = 1; i < segment->parameters.size(); i++) {
            std::cout << "\"" << segment->parameters[i] << "\" ";
        }
        std::cout << std::endl;
    }

    if (!segment->children.empty()) {
        std::cout << "  �ӽڵ�:" << std::endl;
        for (const auto& child : segment->children) {
            std::cout << "    - " << child->name;
            if (!child->parameters.empty()) {
                std::cout << " : ";
                for (const auto& param : child->parameters) {
                    std::cout << "\"" << param << "\" ";
                }
            }
            std::cout << std::endl;
        }
    }
}

// ��ӡvia��Ϣ
void KiCadParser::printViaInfo(const std::shared_ptr<Node>& via) {
    if (!via || via->name != "via") {
        std::cout << "������Ч��via�ڵ�" << std::endl;
        return;
    }

    // �����У���ʾ��ţ��ײ�����
    std::cout << "Via";
    if (!via->parameters.empty()) {
        std::cout << " #" << via->parameters[0];
    }
    std::cout << ":" << std::endl;

    // �����У��ӵڶ���������ʼ��ӡ��������ţ�
    if (via->parameters.size() > 1) {
        std::cout << "  ����: ";
        for (size_t i = 1; i < via->parameters.size(); ++i) {
            std::cout << "\"" << via->parameters[i] << "\" ";
        }
        std::cout << std::endl;
    }

    // �ӽڵ��ӡ���ֲ���
    if (!via->children.empty()) {
        std::cout << "  �ӽڵ�:" << std::endl;
        for (const auto& child : via->children) {
            std::cout << "    - " << child->name;
            if (!child->parameters.empty()) {
                std::cout << " : ";
                for (const auto& param : child->parameters) {
                    std::cout << "\"" << param << "\" ";
                }
            }
            std::cout << std::endl;
        }
    }
}


// ��ӡfootprint��Ϣ
void KiCadParser::printFootprintInfo(const std::shared_ptr<Node>& footprint) {
    if (!footprint || footprint->name != "footprint") {
        std::cout << "������Ч��footprint�ڵ�" << std::endl;
        return;
    }

    std::cout << "Footprint: ";
    if (!footprint->parameters.empty()) {
        std::cout << footprint->parameters[0];
    }
    std::cout << std::endl;

    if (!footprint->children.empty()) {
        std::cout << "  �ӽڵ�:" << std::endl;
        for (const auto& child : footprint->children) {
            std::cout << "    - " << child->name;
            if (!child->parameters.empty()) {
                std::cout << " : ";
                for (const auto& param : child->parameters) {
                    std::cout << "\"" << param << "\" ";
                }
            }
            std::cout << " (�ӽڵ�����: " << child->children.size() << ")" << std::endl;
        }
    }
}

// ��ӡ����segments
void KiCadParser::printAllSegments() {
    if (allSegments.empty()) {
        std::cout << "û���ҵ��κ�segment" << std::endl;
        return;
    }

    std::cout << "\n=== ����Segment�б� (�� " << allSegments.size() << " ��) ===" << std::endl;
    for (const auto& segment : allSegments) {
        printSegmentInfo(segment);
        std::cout << std::endl;
    }
}

// ��ӡ����vias
void KiCadParser::printAllVias() {
    if (allVias.empty()) {
        std::cout << "û���ҵ��κ�via" << std::endl;
        return;
    }

    std::cout << "\n=== ����Via�б� (�� " << allVias.size() << " ��) ===" << std::endl;
    for (const auto& via : allVias) {
        printViaInfo(via);
        std::cout << std::endl;
    }
}

// ��ӡ����footprints
void KiCadParser::printAllFootprints() {
    if (allFootprints.empty()) {
        std::cout << "û���ҵ��κ�footprint" << std::endl;
        return;
    }

    std::cout << "\n=== ����Footprint�б� (�� " << allFootprints.size() << " ��) ===" << std::endl;
    for (const auto& footprint : allFootprints) {
        printFootprintInfo(footprint);
        std::cout << std::endl;
    }
}

// ��ӡ��״�ṹ
void KiCadParser::printStructure(const std::shared_ptr<Node>& node, int depth, int maxDepth) {
    std::shared_ptr<Node> current = node ? node : root;

    if (!current || depth > maxDepth) {
        return;
    }

    std::string indent(depth * 2, ' ');
    std::cout << indent << "Node: " << current->name;

    if (!current->parameters.empty()) {
        std::cout << " | Parameters: ";
        for (const auto& param : current->parameters) {
            std::cout << "\"" << param << "\" ";
        }
    }
    std::cout << std::endl;

    for (const auto& child : current->children) {
        printStructure(child, depth + 1, maxDepth);
    }
}





// ����֮ǰ�� Node �� KiCadParser �඼�Ѱ������ڣ��������µ�������

// �ݹ��ӡ�ڵ�ṹ
void printNodeStructure(const std::shared_ptr<Node>& node, int depth = 0) {
    if (!node) return;



    // �����ǰ�ڵ�����ƺͲ���
    std::string indent(depth * 2, ' ');  // ����
    std::cout << indent << "�ڵ�����: " << node->name << std::endl;

    if (!node->parameters.empty()) {
        std::cout << indent << "  ����: ";
        for (const auto& param : node->parameters) {
            std::cout << "\"" << param << "\" ";
        }
        std::cout << std::endl;
    }

    // �ݹ�����ӽڵ�
    if (!node->children.empty()) {
        std::cout << indent << "  �ӽڵ�:" << std::endl;
        for (const auto& child : node->children) {
            printNodeStructure(child, depth + 1);  // ��ȼ�1���еݹ�
        }
    }
}

// === ������ͨ�õĵݹ��ӡ�����ڴ�ӡ via �������� ===
static void printNodeGeneric(const std::shared_ptr<Node>& node, int depth = 0) {
    if (!node) return;
    std::string indent(depth * 2, ' ');

    // ���ף��ڵ���
    std::cout << indent << "- " << node->name;

    // ͬ�У��ڵ����������У�
    if (!node->parameters.empty()) {
        std::cout << " : ";
        for (const auto& p : node->parameters) {
            std::cout << "\"" << p << "\" ";
        }
    }
    std::cout << std::endl;

    // �ݹ��ӡ�ӽڵ�
    for (const auto& child : node->children) {
        printNodeGeneric(child, depth + 1);
    }
}

// === ר���� via �������ݹ��ӡ���ײ���Ϊ��ţ��� segment/footprint һ�£� ===
static void printViaFull(const std::shared_ptr<Node>& via, int depth = 0) {
    if (!via || via->name != "via") return;
    std::string indent(depth * 2, ' ');

    // ���⣺Via #���
    std::cout << indent << "Via";
    if (!via->parameters.empty()) {
        std::cout << " #" << via->parameters[0]; // Լ�����ײ�Ϊ���
    }
    std::cout << std::endl;

    // ��ӡ���������������������У�
    if (via->parameters.size() > 1) {
        std::cout << indent << "  ����: ";
        for (size_t i = 1; i < via->parameters.size(); ++i) {
            std::cout << "\"" << via->parameters[i] << "\" ";
        }
        std::cout << std::endl;
    }

    // �ݹ��ӡ�������������� at/size/drill/layers/net �ȣ�
    if (!via->children.empty()) {
        std::cout << indent << "  ����:" << std::endl;
        for (const auto& child : via->children) {
            printNodeGeneric(child, depth + 2);
        }
    }
}

// ===================== Edge.Cuts ���γߴ���㣨������ =====================

std::string KiCadParser::unquote(const std::string& s) {
    if (s.size() >= 2 && s.front() == '"' && s.back() == '"') {
        return s.substr(1, s.size() - 2);
    }
    return s;
}

bool KiCadParser::isEdgeCuts(const std::shared_ptr<Node>& n) {
    if (!n) return false;
    for (const auto& ch : n->children) {
        if (ch && ch->name == "layer" && !ch->parameters.empty()) {
            if (unquote(ch->parameters[0]) == "Edge.Cuts") return true;
        }
    }
    return false;
}

bool KiCadParser::readXY(const std::shared_ptr<Node>& n, const char* childName, double& x, double& y) {
    if (!n) return false;
    for (const auto& ch : n->children) {
        if (ch && ch->name == childName) {
            if (ch->parameters.size() >= 2) {
                try {
                    x = std::stod(ch->parameters[0]);
                    y = std::stod(ch->parameters[1]);
                    return true;
                }
                catch (...) {
                    return false;
                }
            }
        }
    }
    return false;
}

inline void KiCadParser::pushPoint(double x, double y, BoardBBox& bbox) {
    if (!bbox.valid) {
        bbox.valid = true;
        bbox.minX = bbox.maxX = x;
        bbox.minY = bbox.maxY = y;
    }
    else {
        if (x < bbox.minX) bbox.minX = x;
        if (x > bbox.maxX) bbox.maxX = x;
        if (y < bbox.minY) bbox.minY = y;
        if (y > bbox.maxY) bbox.maxY = y;
    }
}

KiCadParser::BoardBBox KiCadParser::getBoardBBox() const {
    BoardBBox bbox;

    // ������ȣ�ֻ���Ķ��� gr_*�����κβ㼶�� gr_*�������� footprint ��� fp_*��
    std::function<void(const std::shared_ptr<Node>&)> dfs =
        [&](const std::shared_ptr<Node>& node) {
        if (!node) return;

        const std::string& nm = node->name;

        // ������ gr_*������ fp_*����װ��ĸ�۲����ڰ����Σ�
        const bool isGraphic =
            (nm == "gr_line" || nm == "gr_arc" || nm == "gr_rect" || nm == "gr_poly");

        if (isGraphic && isEdgeCuts(node)) {
            if (nm == "gr_line") {
                double x1, y1, x2, y2;
                if (readXY(node, "start", x1, y1)) pushPoint(x1, y1, bbox);
                if (readXY(node, "end", x2, y2)) pushPoint(x2, y2, bbox);
            }
            else if (nm == "gr_arc") {
                double xs, ys, xe, ye, xm, ym;
                if (readXY(node, "start", xs, ys)) pushPoint(xs, ys, bbox);
                if (readXY(node, "end", xe, ye)) pushPoint(xe, ye, bbox);
                // ��ѡ���� mid Ҳ���룬�ݴ����
                if (readXY(node, "mid", xm, ym)) pushPoint(xm, ym, bbox);
            }
            else if (nm == "gr_rect") {
                double x1, y1, x2, y2;
                if (readXY(node, "start", x1, y1)) pushPoint(x1, y1, bbox);
                if (readXY(node, "end", x2, y2)) pushPoint(x2, y2, bbox);
            }
            else if (nm == "gr_poly") {
                // (gr_poly (pts (xy x y) (xy x y) ...)) ���� ȡ���� xy
                for (const auto& ch : node->children) {
                    if (ch && ch->name == "pts") {
                        for (const auto& pt : ch->children) {
                            if (pt && pt->name == "xy" && pt->parameters.size() >= 2) {
                                try {
                                    double x = std::stod(pt->parameters[0]);
                                    double y = std::stod(pt->parameters[1]);
                                    pushPoint(x, y, bbox);
                                }
                                catch (...) { /* ���Ե������ʧ�� */ }
                            }
                        }
                    }
                }
            }
        }

        // �����ݹ�
        for (const auto& ch : node->children) {
            dfs(ch);
        }
        };

    // �Ӹ��ڵ���������������� getAllNodes ����ڵ�������Ϊ���׵ݹ���������
    // ������� root����ֱ�� dfs(root)���������������нӿڣ�
    const auto& allNodes = getAllNodes();
    for (const auto& n : allNodes) dfs(n);

    return bbox;
}

bool KiCadParser::getBoardSize(double& width, double& height) const {
    auto bb = getBoardBBox();
    if (!bb.valid) return false;
    width = bb.width();
    height = bb.height();
    return true;
}

bool KiCadParser::getBoardCorners(Point2D& bl,
    Point2D& tl,
    Point2D& tr,
    Point2D& br) const
{
    auto bb = getBoardBBox();            // ������ʵ�ֵ���Ӿ���
    if (!bb.valid) return false;

    bl = { bb.minX, bb.minY };           // ����
    tl = { bb.minX, bb.maxY };           // ����
    tr = { bb.maxX, bb.maxY };           // ����
    br = { bb.maxX, bb.minY };           // ����
    return true;
}


#include <iostream>
#include <iomanip>
#include <string>
#include "First_Part12.h"
/*
int main(int argc, char** argv) {
    auto runOnce = [](const std::string& path) {
        KiCadParser parser;
        if (!parser.parseFile(path)) {
            std::cerr << "[����ʧ��] " << path << std::endl;
            return;
        }

        // ����Ľ�
        KiCadParser::Point2D bl, tl, tr, br;
        if (!parser.getBoardCorners(bl, tl, tr, br)) {
            std::cout << "[δ�ҵ� Edge.Cuts ���] " << path << std::endl;
            return;
        }

        // ��ߣ����� bbox ���㣬Ҳ��ֱ���� getBoardSize��
        double w = 0.0, h = 0.0;
        parser.getBoardSize(w, h);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "�ļ�: " << path << "\n";
        std::cout << "����Ľǵ㣨��λ: mm��\n";
        std::cout << "  BL(����): (" << bl.x << ", " << bl.y << ")\n";
        std::cout << "  TL(����): (" << tl.x << ", " << tl.y << ")\n";
        std::cout << "  TR(����): (" << tr.x << ", " << tr.y << ")\n";
        std::cout << "  BR(����): (" << br.x << ", " << br.y << ")\n";
        std::cout << "�ߴ磺��� = " << w << "���߶� = " << h << "\n";
        std::cout << "---------------------------------------------\n";
        };

    if (argc >= 2) {
        for (int i = 1; i < argc; ++i) runOnce(argv[i]);
    }
    else {
        // �ĳ������������·�������� L92.txt / A30.txt / B40.txt
        runOnce(R"(D:\����ʫ����Ŀ\Liang�Ŀ���\7.0.0�汾\L92.txt)");
    }

    return 0;
}
*/