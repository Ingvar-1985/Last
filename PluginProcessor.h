#pragma once

#include <JuceHeader.h>
#include "JCM2203Preamp.h"
#include <array>
#include <memory>

class JCM2203ModdedPreampAudioProcessor : public juce::AudioProcessor
{
public:
    JCM2203ModdedPreampAudioProcessor();
    ~JCM2203ModdedPreampAudioProcessor() override = default;

    void prepareToPlay(double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

#ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported(const BusesLayout& layouts) const override;
#endif

    void processBlock(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
    void processBlock(juce::AudioBuffer<double>&, juce::MidiBuffer&) override;

    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override { return true; }

    const juce::String getName() const override { return JucePlugin_Name; }
    bool acceptsMidi() const override { return false; }
    bool producesMidi() const override { return false; }
    bool isMidiEffect() const override { return false; }
    double getTailLengthSeconds() const override { return 0.1; }

    int getNumPrograms() override { return 1; }
    int getCurrentProgram() override { return 0; }
    void setCurrentProgram(int) override {}
    const juce::String getProgramName(int) override { return "Default"; }
    void changeProgramName(int, const juce::String&) override {}

    void getStateInformation(juce::MemoryBlock& destData) override;
    void setStateInformation(const void* data, int sizeInBytes) override;

    juce::AudioProcessorValueTreeState apvts;

private:
    static juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout();

    template <typename SampleType>
    void processInternal(juce::AudioBuffer<SampleType>& buffer);

    void updateDspParameters(bool immediate = false);

    // 8x Oversampling (3 stages of 2x FIR equiripple in double precision)
    std::unique_ptr<juce::dsp::Oversampling<double>> oversampling;

    // One independent physical tube preamp circuit solver per channel (max stereo)
    std::array<jcm2203::JCM2203PreampChannel, 2> channels;
    juce::AudioBuffer<double> doubleScratchBuffer;

    // Host-rate parameter smoothers (combined with internal oversampled smoothing)
    juce::SmoothedValue<double, juce::ValueSmoothingTypes::Linear> gain1Smoothed;
    juce::SmoothedValue<double, juce::ValueSmoothingTypes::Linear> gain2Smoothed;
    juce::SmoothedValue<double, juce::ValueSmoothingTypes::Linear> trebleSmoothed;
    juce::SmoothedValue<double, juce::ValueSmoothingTypes::Linear> bassSmoothed;
    juce::SmoothedValue<double, juce::ValueSmoothingTypes::Linear> middleSmoothed;
    juce::SmoothedValue<double, juce::ValueSmoothingTypes::Linear> masterSmoothed;
    juce::SmoothedValue<double, juce::ValueSmoothingTypes::Linear> inputVoltsSmoothed;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(JCM2203ModdedPreampAudioProcessor)
};
