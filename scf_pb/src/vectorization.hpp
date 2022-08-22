#pragma once

#include<tuple>
#include<vector>

#define MAKE_VECTORIZED_FUNCTION(f) auto f##_v = cartesian_product::vectorize_function(f)

namespace cartesian_product{
// cross_imp(f, v...) means "do `f` for each element of cartesian product of v..."
template<typename F>
inline void cross_imp(F f) {
    f();
}

template<typename F, typename H, typename... Ts>
inline void cross_imp(F f, std::vector<H> const& h,
                           std::vector<Ts> const&... t) {
    for(H const& he: h)
        cross_imp([&](Ts const&... ts){
                      f(he, ts...);
                  }, t...);
}

template<typename... Ts>
std::vector<std::tuple<Ts...>> cross(std::vector<Ts> const&... in) {
    std::vector<std::tuple<Ts...>> res;
    cross_imp([&](Ts const&... ts){
                  res.emplace_back(ts...);
              }, in...);
    return res;
}


template<typename F, typename... Ts>
auto cross_apply(F f, std::vector<Ts> const&... in) {
    typedef decltype(f(in[0]...)) ReturnType;
    std::vector<ReturnType> res;
    cross_imp([&](Ts const&... ts){
                  res.emplace_back(f(ts...));
              }, in...);
    return res;
}

template <typename F>
auto vectorize_function(F f){
    auto vectorized_func = [f](auto const&... args){return cross_apply(f, args...);};
    return vectorized_func;
}

}
