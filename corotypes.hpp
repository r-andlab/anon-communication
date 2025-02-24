#ifndef __COROTYPES_HPP__
#define __COROTYPES_HPP__

#include <functional>
#include <boost/coroutine2/coroutine.hpp>

using coro_t = boost::coroutines2::coroutine<void>::pull_type;
using yield_t = boost::coroutines2::coroutine<void>::push_type;
using coro_lambda_t = std::function<void(yield_t&)>;

#endif
