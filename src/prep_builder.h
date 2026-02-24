#pragma once

#include <Rcpp.h>
#include <memory>

#include "context.h"
#include "proto.h"

namespace uuber {

NativePrepProto build_prep_proto(const Rcpp::List& prep);
std::unique_ptr<NativeContext> build_context_from_proto(const NativePrepProto& proto);

}
