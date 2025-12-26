#pragma once

/* This file was generated using schema-salad code generator.
 *
 * The embedded document is subject to the license of the original schema.
 */

#include <any>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace example_com {

struct store_config {
    bool simplifyTypes = true;
    bool transformListsToMaps = true;
    bool generateTags = false;
};

inline auto simplifyType(YAML::Node type, store_config const& config) -> YAML::Node {
    if (!config.simplifyTypes) return type;
    auto is_optional = [](YAML::Node const & node) {
        return node.IsSequence() && node.size() == 2u && node[0].Scalar() == "null";
    };

    auto is_array = [](YAML::Node const & node) {
        return node.IsMap() && node["type"].Scalar() == "array" && node["items"].IsScalar();
    };

    // 1. Collapsing optional scalar types into one option
    if (is_optional(type) && type[1].IsScalar()) {
        type = type[1].as<std::string>() + "?";
    }

    // 2. Collapsing array types into one option
    if (is_array(type)) {
        type = type["items"].as<std::string>() + "[]";
    }

    // 3. Collapsing optional array types into one option
    if (is_optional(type) && is_array(type[1])) {
        type = type[1]["items"].as<std::string>() + "[]?";
    }

    return type;
}

inline auto expandType(YAML::Node type) -> YAML::Node {
    auto ends_with = [](std::string str, std::string suffix) {
        if (str.size() < suffix.size()) return false;
        auto str_suffix = str.substr(str.size()-suffix.size(), suffix.size());
        return str_suffix == suffix;
    };

    // 0. If not a scalar type, nothing to do
    if (!type.IsDefined() || !type.IsScalar()) {
        return type;
    }

    auto str = type.as<std::string>();
    // 1. Check if optional array type and expand
    if (ends_with(str, "[]?")) {
        auto result = YAML::Node{};
        result.push_back(YAML::Node{"null"});
        auto array = YAML::Node{};
        array["type"] = "array";
        array["items"] = expandType(YAML::Node(str.substr(0, str.size()-3)));
        result.push_back(array);
        return result;
    }

    // 2. Expand array
    if (ends_with(str, "[]")) {
        auto array = YAML::Node{};
        array["type"] = "array";
        array["items"] = expandType(YAML::Node(str.substr(0, str.size()-2)));
        return array;
    }

    // 3. Expand optional scalar type
    if (ends_with(str, "?")) {
        auto result = YAML::Node{};
        result.push_back(YAML::Node{"null"});
        result.push_back(expandType(YAML::Node(str.substr(0, str.size()-1))));
        return result;
    }
    return type;
}

inline auto mergeYaml(YAML::Node n1, YAML::Node n2) {
    for (auto const& e : n2) {
        n1[e.first.as<std::string>()] = e.second;
    }
    return n1;
}

// declaring toYaml
inline auto toYaml(bool v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(float v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(double v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(char v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(int8_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(uint8_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(int16_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(uint16_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(int32_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(uint32_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(int64_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(uint64_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(std::monostate const&, [[maybe_unused]] store_config const&) {
    return YAML::Node(YAML::NodeType::Undefined);
}
inline auto toYaml(std::string const& v, [[maybe_unused]] store_config const&) {
    return YAML::Node{v};
}

template <typename T, typename ...Args>
auto anyToYaml_impl(std::any const& a, [[maybe_unused]] store_config const& config) {
    if (auto v = std::any_cast<T const>(&a)) {
        return toYaml(*v, config);
    }
    if constexpr (sizeof...(Args) > 0) {
        return anyToYaml_impl<Args...>(a, config);
    }
    return toYaml(std::monostate{}, config);
}

inline auto toYaml(std::any const& a, [[maybe_unused]] store_config const& config) {
    return anyToYaml_impl<bool,
                          float,
                          double,
                          char,
                          int8_t,
                          uint8_t,
                          int16_t,
                          uint16_t,
                          int32_t,
                          uint32_t,
                          int64_t,
                          uint64_t,
                          std::string>(a, config);
}

// declaring fromYaml
inline void fromYaml(YAML::Node const& n, bool& v) {
    v = n.as<bool>();
}
inline void fromYaml(YAML::Node const& n, float& v) {
    v = n.as<float>();
}
inline void fromYaml(YAML::Node const& n, double& v) {
    v = n.as<double>();
}
inline void fromYaml(YAML::Node const& n, int32_t& v) {
    v = n.as<int32_t>();
}
inline void fromYaml(YAML::Node const& n, int64_t& v) {
    v = n.as<int64_t>();
}
inline void fromYaml(YAML::Node const& n, std::string& v) {
    v = n.as<std::string>();
}
inline void fromYaml(YAML::Node const&, std::any&) {
}
inline void fromYaml(YAML::Node const&, std::monostate&) {
}

inline void addYamlField(YAML::Node& node, std::string const& key, YAML::Node value) {
    if (value.IsDefined()) {
        node[key] = value;
    }
}

inline auto convertListToMap(YAML::Node list, std::string const& mapSubject,
                             std::string const& mapPredicate, store_config const& config) {
    if (!config.transformListsToMaps) return list;
    if (mapSubject.empty()) return list;
    if (list.size() == 0) return list;
    auto map = YAML::Node{};
    for (YAML::Node n : list) {
        auto key = n[mapSubject].as<std::string>();
        if (mapPredicate.empty() || n[mapPredicate].IsMap() || n.size() > 2) {
            n.remove(mapSubject);
            map[key] = n;
        } else {
            map[key] = n[mapPredicate];
        }
    }
    return map;
}
inline auto convertMapToList(YAML::Node map, std::string const& mapSubject,
                             std::string const& mapPredicate) {
    if (mapSubject.empty()) return map;
    if (!map.IsDefined()) return map;
    if (!map.IsMap()) return map;
    auto list = YAML::Node{};
    for (auto n : map) {
        if (mapPredicate.empty() || n.second.IsMap()) {
            n.second[mapSubject] = n.first;
            list.push_back(n.second);
        } else {
            auto n2 = YAML::Node{};
            n2[mapSubject] = n.first;
            n2[mapPredicate] = n.second;
            list.push_back(n2);
        }
    }
    return list;
}

template <typename T> struct IsConstant : std::false_type {};

// fwd declaring toYaml
template <typename T>
auto toYaml(std::vector<T> const& v, [[maybe_unused]] store_config const& config) -> YAML::Node;
template <typename T>
auto toYaml(std::map<std::string, T> const& v, [[maybe_unused]] store_config const& config) -> YAML::Node;
template <typename T>
auto toYaml(T const& t, [[maybe_unused]] store_config const& config) -> YAML::Node;
template <typename ...Args>
auto toYaml(std::variant<Args...> const& t, [[maybe_unused]] store_config const& config) -> YAML::Node;

// fwd declaring fromYaml
template <typename T>
void fromYaml(YAML::Node const& n, std::vector<T>& v);
template <typename T>
void fromYaml(YAML::Node const& n, std::map<std::string, T>& v);
template <typename T>
void fromYaml(YAML::Node const& n, T& t);
template <typename ...Args>
void fromYaml(YAML::Node const& n, std::variant<Args...>& t);

template <typename T>
struct DetectAndExtractFromYaml {
    auto operator()(YAML::Node const&) const -> std::optional<T> {
        return std::nullopt;
    }
};

// special cwl expression string
struct cwl_expression_string {
    std::string s;

    auto toYaml([[maybe_unused]] store_config const& config) const {
        auto n = YAML::Node{s};
        if (config.generateTags) {
            n.SetTag("Expression");
        }
        return n;
    }
    void fromYaml(YAML::Node const& n) {
        s = n.as<std::string>();
    }
};


template <>
struct DetectAndExtractFromYaml<std::monostate> {
    auto operator()(YAML::Node const& n) const -> std::optional<std::monostate> {
        if (!n.IsDefined()) return std::monostate{};
        return std::nullopt;
    }
};

template <typename S>
struct DetectAndExtractFromYaml_implScalar {
    auto operator()(YAML::Node const& n) const -> std::optional<S> {
        try {
            if (n.IsScalar()) return n.as<S>();
        } catch(...) {}
        return std::nullopt;
    }
};

template <> struct DetectAndExtractFromYaml<bool>        : DetectAndExtractFromYaml_implScalar<bool>{};
template <> struct DetectAndExtractFromYaml<float>       : DetectAndExtractFromYaml_implScalar<float>{};
template <> struct DetectAndExtractFromYaml<double>      : DetectAndExtractFromYaml_implScalar<double>{};
template <> struct DetectAndExtractFromYaml<int32_t>     : DetectAndExtractFromYaml_implScalar<int32_t>{};
template <> struct DetectAndExtractFromYaml<int64_t>     : DetectAndExtractFromYaml_implScalar<int64_t>{};
template <> struct DetectAndExtractFromYaml<std::string> : DetectAndExtractFromYaml_implScalar<std::string>{};

template <typename T>
struct DetectAndExtractFromYaml<std::vector<T>> {
    auto operator()(YAML::Node const& n) const -> std::optional<std::vector<T>> {
        if (!n.IsDefined()) return std::nullopt;
        if (!n.IsSequence()) return std::nullopt;
        auto res = std::vector<T>{};
        fromYaml(n, res);
        return res;
    }
};

template <typename T>
struct DetectAndExtractFromYaml<std::map<std::string, T>> {
    auto operator()(YAML::Node const& n) const -> std::optional<std::map<std::string, T>> {
        if (!n.IsDefined()) return std::nullopt;
        if (!n.IsMap()) return std::nullopt;
        auto res = std::map<std::string, T>{};
        fromYaml(n, res);
        return res;
    }
};

template <typename T>
class heap_object {
    std::unique_ptr<T> data = std::make_unique<T>();

public:
    using value_t = T;
    heap_object() noexcept(false) = default;
    heap_object(heap_object const& oth) {
        *data = *oth;
    }
    heap_object(heap_object&& oth) noexcept(noexcept(*data = std::move(*oth))) {
        *data = std::move(*oth);
    }

    template <typename T2>
    heap_object(T2 const& oth) {
        *data = oth;
    }
    template <typename T2>
    heap_object(T2&& oth) noexcept(noexcept(*data = std::forward<T2>(oth))) {
        *data = std::forward<T2>(oth);
    }

    ~heap_object();

    auto operator=(heap_object const& oth) -> heap_object& {
        *data = *oth;
        return *this;
    }
    auto operator=(heap_object&& oth) noexcept(noexcept(*data = std::move(*oth))) -> heap_object& {
        *data = std::move(*oth);
        return *this;
    }

    template <typename T2>
    auto operator=(T2 const& oth) -> heap_object& {
        *data = oth;
        return *this;
    }
    template <typename T2>
    auto operator=(T2&& oth) noexcept(noexcept(*data = std::forward<T2>(oth))) -> heap_object& {
        *data = std::forward<T2>(oth);
        return *this;
    }

    auto operator->() noexcept(true) -> T* {
        return data.get();
    }
    auto operator->() const noexcept(true) -> T const* {
        return data.get();
    }
    auto operator*() noexcept(true) -> T& {
        return *data;
    }
    auto operator*() const noexcept(true) -> T const& {
        return *data;
    }
};

}
namespace example_com { struct MyRecordOne; }
namespace example_com { struct MyRecordTwo; }
namespace example_com {
struct MyRecordOne {
    heap_object<std::string> name;
    virtual ~MyRecordOne() = default;
    virtual auto toYaml([[maybe_unused]] example_com::store_config const& config) const -> YAML::Node;
    virtual void fromYaml(YAML::Node const& n);
};
}

namespace example_com {
struct MyRecordTwo
    : example_com::MyRecordOne {
    heap_object<int32_t> value;
    ~MyRecordTwo() override = default;
    auto toYaml([[maybe_unused]] example_com::store_config const& config) const -> YAML::Node override;
    void fromYaml(YAML::Node const& n) override;
};
}

namespace example_com {
template <typename T> heap_object<T>::~heap_object() = default;
}

inline auto example_com::MyRecordOne::toYaml([[maybe_unused]] ::example_com::store_config const& config) const -> YAML::Node {
    using ::example_com::toYaml;
    auto n = YAML::Node{};
    if (config.generateTags) {
        n.SetTag("MyRecordOne");
    }
    {
         auto member = toYaml(*name, config);
         member = convertListToMap(member, "", "", config);
        addYamlField(n, "name", member);
    }
    return n;
}
inline void example_com::MyRecordOne::fromYaml([[maybe_unused]] YAML::Node const& n) {
    using ::example_com::fromYaml;
    {
        auto nodeAsList = convertMapToList(n["name"], "", "");
        auto expandedNode = (nodeAsList);
        fromYaml(expandedNode, *name);
    }
}
namespace example_com {
template <>
struct DetectAndExtractFromYaml<::example_com::MyRecordOne> {
    auto operator()(YAML::Node const& n) const -> std::optional<::example_com::MyRecordOne> {
        if (!n.IsDefined()) return std::nullopt;
        if (!n.IsMap()) return std::nullopt;
        auto res = ::example_com::MyRecordOne{};

        if constexpr (::example_com::IsConstant<decltype(res.name)::value_t>::value) try {
            fromYaml(n["name"], *res.name);
            fromYaml(n, res);
            return res;
        } catch(...) {}

        return std::nullopt;
    }
};
}
inline auto example_com::MyRecordTwo::toYaml([[maybe_unused]] ::example_com::store_config const& config) const -> YAML::Node {
    using ::example_com::toYaml;
    auto n = YAML::Node{};
    if (config.generateTags) {
        n.SetTag("MyRecordTwo");
    }
    n = mergeYaml(n, example_com::MyRecordOne::toYaml(config));
    {
         auto member = toYaml(*value, config);
         member = convertListToMap(member, "", "", config);
        addYamlField(n, "value", member);
    }
    return n;
}
inline void example_com::MyRecordTwo::fromYaml([[maybe_unused]] YAML::Node const& n) {
    using ::example_com::fromYaml;
    example_com::MyRecordOne::fromYaml(n);
    {
        auto nodeAsList = convertMapToList(n["value"], "", "");
        auto expandedNode = (nodeAsList);
        fromYaml(expandedNode, *value);
    }
}
namespace example_com {
template <>
struct DetectAndExtractFromYaml<::example_com::MyRecordTwo> {
    auto operator()(YAML::Node const& n) const -> std::optional<::example_com::MyRecordTwo> {
        if (!n.IsDefined()) return std::nullopt;
        if (!n.IsMap()) return std::nullopt;
        auto res = ::example_com::MyRecordTwo{};

        if constexpr (::example_com::IsConstant<decltype(res.value)::value_t>::value) try {
            fromYaml(n["value"], *res.value);
            fromYaml(n, res);
            return res;
        } catch(...) {}

        return std::nullopt;
    }
};
}
namespace example_com {

template <typename T>
auto toYaml(std::vector<T> const& v, [[maybe_unused]] store_config const& config) -> YAML::Node {
    auto n = YAML::Node(YAML::NodeType::Sequence);
    for (auto const& e : v) {
        n.push_back(toYaml(e, config));
    }
    return n;
}

template <typename T>
auto toYaml(std::map<std::string, T> const& v, [[maybe_unused]] store_config const& config) -> YAML::Node {
    auto n = YAML::Node(YAML::NodeType::Map);
    for (auto const& [key, value] : v) {
        n[key] = toYaml(value, config);
    }
    return n;
}

template <typename T>
auto toYaml(T const& t, [[maybe_unused]] store_config const& config) -> YAML::Node {
    if constexpr (std::is_enum_v<T>) {
        return toYaml(t, config);
    } else {
        return t.toYaml(config);
    }
}

template <typename ...Args>
auto toYaml(std::variant<Args...> const& t, store_config const& config) -> YAML::Node {
    return std::visit([config](auto const& e) {
        return toYaml(e, config);
    }, t);
}

template <typename T>
void fromYaml(YAML::Node const& n, std::vector<T>& v){
    if (!n.IsSequence()) return;
    for (auto e : n) {
        v.emplace_back();
        fromYaml(e, v.back());
    }
}

template <typename T>
void fromYaml(YAML::Node const& n, std::map<std::string, T>& v){
    if (!n.IsMap()) return;
    for (auto e : n) {
        auto key = e.first.as<std::string>();
        fromYaml(e.second, v[key]);
    }
}

template <typename T>
void fromYaml(YAML::Node const& n, T& t){
    if constexpr (std::is_enum_v<T>) {
        fromYaml(n, t);
    } else {
        t.fromYaml(n);
    }
}

template <typename SomeVariant, typename Head, typename ...Args>
bool detectAndExtractFromYaml(YAML::Node const& n, SomeVariant& v, Head* = nullptr) {
    auto r = DetectAndExtractFromYaml<Head>{}(n);
    if (r) {
        v = *r;
        return true;
    }
    if constexpr (sizeof...(Args) > 0) {
        return detectAndExtractFromYaml<SomeVariant, Args...>(n, v);
    }
    return false;
}

template <typename SomeVariant, typename Head, typename Tail>
bool detectAndExtractFromYaml(YAML::Node const& n, std::variant<std::monostate, Tail>& v, Head* = nullptr) {
    auto r = DetectAndExtractFromYaml<Head>{}(n);
    if (r) {
        v = *r;
        return true;
    }
    auto t = Tail{};
    fromYaml(n, t);
    v = t;
    return true;
}

template <typename ...Args>
void fromYaml(YAML::Node const& n, std::variant<Args...>& v){
    bool found = detectAndExtractFromYaml<std::variant<Args...>, Args...>(n, v);
    if (!found) throw std::runtime_error{"didn't find any overload"};
}
using DocumentRootType = std::variant<>;
auto load_document_from_yaml(YAML::Node n) -> DocumentRootType {
    DocumentRootType root;
    fromYaml(n, root);
    return root;
}
auto load_document_from_string(std::string document) -> DocumentRootType {
    return load_document_from_yaml(YAML::Load(document));
}
auto load_document(std::filesystem::path path) -> DocumentRootType {
    return load_document_from_yaml(YAML::LoadFile(path.string()));
}
void store_document(DocumentRootType const& root, std::ostream& ostream, store_config config={}) {
    auto y = toYaml(root, config);

    YAML::Emitter out;
    out << y;
    ostream << out.c_str() << std::endl;
}
void store_document(DocumentRootType const& root, std::filesystem::path const& path, store_config config={}) {
    auto ofs = std::ofstream{path};
    store_document(root, ofs, config);
}
auto store_document_as_string(DocumentRootType const& root, store_config config={}) -> std::string {
    auto ss = std::stringstream{};
    store_document(root, ss, config);
    return ss.str();
}

}