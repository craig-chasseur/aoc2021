load("@rules_cc//cc:defs.bzl", "cc_binary")

cc_binary(
    name = "part1",
    srcs = ["part1.cc"],
    deps = [
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_googlesource_code_re2//:re2",
        "//util:check",
        "//util:grid2",
        "//util:io",
    ],
)
