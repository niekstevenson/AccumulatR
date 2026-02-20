#pragma once

#include "context.h"

namespace uuber {

// Compile dense IR nodes to a single topologically ordered opcode program.
KernelProgram compile_kernel_program(const IrContext &ir);

} // namespace uuber

