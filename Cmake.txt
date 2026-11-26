cmake_minimum_required(VERSION 3.22)
project(JCM2203ModdedPreamp VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Універсальна збірка для macOS (Intel + Apple Silicon)
set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE INTERNAL "")

include(FetchContent)
FetchContent_Declare(
    JUCE
    GIT_REPOSITORY https://github.com
    GIT_TAG        8.0.4
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(JUCE)

juce_add_plugin(JCM2203ModdedPreamp
    COMPANY_NAME "Valkyrie DSP"
    BUNDLE_ID "com.valkyriedsp.jcm2203preamp"
    IS_SYNTH FALSE
    NEEDS_MIDI_INPUT FALSE
    NEEDS_MIDI_OUTPUT FALSE
    IS_MIDI_EFFECT FALSE
    EDITOR_WANTS_KEYBOARD_FOCUS FALSE
    COPY_PLUGIN_AFTER_BUILD FALSE
    PLUGIN_MANUFACTURER_CODE Vdsp
    PLUGIN_CODE Jc23
    FORMATS VST3
    PRODUCT_NAME "JCM2203 Modded Preamp"
)

juce_generate_juce_header(JCM2203ModdedPreamp)

target_sources(JCM2203ModdedPreamp
    PRIVATE
        JCM2203Preamp.h
        JCM2203Preamp.cpp
        PluginProcessor.h
        PluginProcessor.cpp
)

target_compile_definitions(JCM2203ModdedPreamp
    PUBLIC
        JUCE_WEB_BROWSER=0
        JUCE_USE_CURL=0
        JUCE_VST3_CAN_REPLACE_VST2=0
)

target_link_libraries(JCM2203ModdedPreamp
    PRIVATE
        juce::juce_audio_utils
        juce::juce_dsp
    PUBLIC
        juce::juce_recommended_config_flags
        juce::juce_recommended_lto_flags
        juce::juce_recommended_warning_flags
)
