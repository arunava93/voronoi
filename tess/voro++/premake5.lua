workspace "VoroGen"

  location "Generated"

  language "C++"
  cppdialect "c++17"

  architecture "x86_64"

  configurations { "Debug", "Release", "Windows" }

  filter { "configurations:Debug" }
    symbols "On"
    linkoptions { "-static" }

  filter { "configurations:Release" }
    optimize "On"
    linkoptions { "-static" }

  filter { "configurations:Windows" }
    optimize "On"
    linkoptions { "-static" }
    makesettings { "CC=x86_64-w64-mingw32-gcc" }
    makesettings { "CXX=x86_64-w64-mingw32-g++" }
    makesettings { "AR=x86_64-w64-mingw32-ar" }

  filter {}

  targetdir ("Build/Bin/%{prj.name}/%{cfg.longname}")
  objdir ("Build/Obj/%{prj.name}/%{cfg.longname}")

function includeEigen()

  includedirs "Libraries/Eigen/"

end

function includeCxxopts()

  includedirs "Libraries/cxxopts/"

end

function useVoroLib()

  includedirs "Projects/voro++/"
  links "VoroLib"

end

project "VoroLib"

  kind "StaticLib"

  files "Projects/voro++/**"
  removefiles {"Projects/voro++/v_base_wl.cc","Projects/voro++/voro++.cc","Projects/voro++/cmd_line.cc"}

project "VoroExe"

  kind "ConsoleApp"

  files {"Projects/voro++/voro++.cc","Projects/voro++/cmd_line.cc"}

  useVoroLib()

project "VoroGen"

  kind "ConsoleApp"
  files "Projects/VoroGen/**"

  includedirs "Projects/voro++"
  includeCxxopts()

  useVoroLib()