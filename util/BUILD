load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "check",
    hdrs = ["check.h"],
)

cc_library(
    name = "dimension_grid",
    hdrs = ["dimension_grid.h"],
    deps = [
        ":check",
    ],
)

cc_library(
    name = "grid2",
    srcs = ["grid2.cc"],
    hdrs = ["grid2.h"],
    deps = [
        ":check",
        "@com_google_absl//absl/strings",
    ],
)

cc_library(
    name = "io",
    srcs = ["io.cc"],
    hdrs = ["io.h"],
    deps = [
        ":check",
        "@com_google_absl//absl/strings",
    ],
)
