load("@rules_cc//cc:defs.bzl", "cc_binary")

cc_binary(
    name = "part1",
    srcs = ["part1.cc"],
    deps = [
        "@com_googlesource_code_re2//:re2",
        "//util:check",
        "//util:io",
    ],
)

cc_binary(
    name = "part2",
    srcs = ["part2.cc"],
    deps = [
        "@com_googlesource_code_re2//:re2",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
        "//util:check",
        "//util:io",
    ],
)
