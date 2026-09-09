{
  "targets": [
    {
      "target_name": "forces",
      "sources": [ "forces.cpp" ],
      "include_dirs": [
        "<(module_root_dir)",
        "<!@(node -p \"require('node-addon-api').include\")"
      ],
      "dependencies": [
        "<!(node -p \"require('node-addon-api').gyp\")"
      ],
      "library_dirs": [
        "<(module_root_dir)"
      ],
      "libraries": [
        "OpenCL.lib"
      ],
      "cflags!": [ "-fno-exceptions" ],
      "cflags_cc!": [ "-fno-exceptions" ],
      "defines": [ "NAPI_DISABLE_CPP_EXCEPTIONS" ]
    }
  ]
}