#include "PluginProcessor.h"

namespace
{
constexpr int kOversamplingFactorOrder = 3; // 2^3 = 8x Oversampling
}

JCM2203ModdedPreampAudioProcessor::JCM2203ModdedPreampAudioProcessor()
    : AudioProcessor(BusesProperties()
                         .withInput("Input", juce::AudioChannelSet::stereo(), true)
                         .withOutput("Output", juce::AudioChannelSet::stereo(), true)),
      apvts(*this, nullptr, "PARAMETERS", createParameterLayout())
{
}

juce::AudioProcessorValueTreeState::ParameterLayout
JCM2203ModdedPreampAudioProcessor::createParameterLayout()
{
    std::vector<std::unique_ptr<juce::RangedAudioParameter>> params;

    // Defaults calibrated to exact XML controlPosition values
    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        juce::ParameterID("gain1", 1), "GAIN 1",
        juce::NormalisableRange<float>(0.0f, 1.0f, 0.001f), 0.896f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        juce::ParameterID("gain2", 1), "GAIN 2",
        juce::NormalisableRange<float>(0.0f, 1.0f, 0.001f), 0.760f));

    params.push_back(std::make_unique<juce::AudioParameterBool>(
        juce::ParameterID("high1", 1), "High 1 (2.2nF)", false));

    params.push_back(std::make_unique<juce::AudioParameterBool>(
        juce::ParameterID("high2", 1), "High 2 (2.2pF)", false));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        juce::ParameterID("treble", 1), "TREBLE",
        juce::NormalisableRange<float>(0.0f, 1.0f, 0.001f), 0.460f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        juce::ParameterID("bass", 1), "BASS",
        juce::NormalisableRange<float>(0.0f, 1.0f, 0.001f), 0.608f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        juce::ParameterID("middle", 1), "MIDDLE",
        juce::NormalisableRange<float>(0.0f, 1.0f, 0.001f), 0.500f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        juce::ParameterID("master", 1), "MASTER",
        juce::NormalisableRange<float>(0.0f, 1.0f, 0.001f), 0.420f));

    // Calibration control: 0 dBFS peak -> Volts (0.15V single coil .. 1.50V active/boosted)
    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        juce::ParameterID("inputVolts", 1), "Input Trim (V @ 0dBFS)",
        juce::NormalisableRange<float>(0.10f, 2.50f, 0.01f, 0.6f), 0.80f));

    return { params.begin(), params.end() };
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool JCM2203ModdedPreampAudioProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
    const auto mainOut = layouts.getMainOutputChannelSet();
    if (mainOut != juce::AudioChannelSet::mono() && mainOut != juce::AudioChannelSet::stereo())
        return false;
    return layouts.getMainInputChannelSet() == mainOut;
}
#endif

void JCM2203ModdedPreampAudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock)
{
    const auto totalOutChannels = static_cast<size_t>(std::max(1, getTotalNumOutputChannels()));

    oversampling = std::make_unique<juce::dsp::Oversampling<double>>(
        totalOutChannels,
        kOversamplingFactorOrder,
        juce::dsp::Oversampling<double>::filterHalfBandFIREquiripple,
        true);
    oversampling->initProcessing(static_cast<size_t>(samplesPerBlock));
    setLatencySamples(static_cast<int>(oversampling->getLatencyInSamples()));

    doubleScratchBuffer.setSize(static_cast<int>(totalOutChannels), samplesPerBlock, false, true, true);

    const double osRate = sampleRate * static_cast<double>(oversampling->getOversamplingFactor());

    gain1Smoothed.reset(sampleRate, 0.015);
    gain2Smoothed.reset(sampleRate, 0.015);
    trebleSmoothed.reset(sampleRate, 0.015);
    bassSmoothed.reset(sampleRate, 0.015);
    middleSmoothed.reset(sampleRate, 0.015);
    masterSmoothed.reset(sampleRate, 0.015);
    inputVoltsSmoothed.reset(sampleRate, 0.015);

    updateDspParameters(true);

    for (auto& ch : channels)
        ch.prepare(osRate);
}

void JCM2203ModdedPreampAudioProcessor::releaseResources()
{
    if (oversampling != nullptr)
        oversampling->reset();
}

void JCM2203ModdedPreampAudioProcessor::updateDspParameters(bool immediate)
{
    const double g1 = apvts.getRawParameterValue("gain1")->load();
    const double g2 = apvts.getRawParameterValue("gain2")->load();
    const double tr = apvts.getRawParameterValue("treble")->load();
    const double ba = apvts.getRawParameterValue("bass")->load();
    const double mi = apvts.getRawParameterValue("middle")->load();
    const double ma = apvts.getRawParameterValue("master")->load();
    const double iv = apvts.getRawParameterValue("inputVolts")->load();

    if (immediate)
    {
        gain1Smoothed.setCurrentAndTargetValue(g1);
        gain2Smoothed.setCurrentAndTargetValue(g2);
        trebleSmoothed.setCurrentAndTargetValue(tr);
        bassSmoothed.setCurrentAndTargetValue(ba);
        middleSmoothed.setCurrentAndTargetValue(mi);
        masterSmoothed.setCurrentAndTargetValue(ma);
        inputVoltsSmoothed.setCurrentAndTargetValue(iv);
    }
    else
    {
        gain1Smoothed.setTargetValue(g1);
        gain2Smoothed.setTargetValue(g2);
        trebleSmoothed.setTargetValue(tr);
        bassSmoothed.setTargetValue(ba);
        middleSmoothed.setTargetValue(mi);
        masterSmoothed.setTargetValue(ma);
        inputVoltsSmoothed.setTargetValue(iv);
    }
}

void JCM2203ModdedPreampAudioProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    processInternal(buffer);
}

void JCM2203ModdedPreampAudioProcessor::processBlock(juce::AudioBuffer<double>& buffer, juce::MidiBuffer&)
{
    processInternal(buffer);
}

template <typename SampleType>
void JCM2203ModdedPreampAudioProcessor::processInternal(juce::AudioBuffer<SampleType>& buffer)
{
    juce::ScopedNoDenormals noDenormals;

    const int numSamples = buffer.getNumSamples();
    const int numIn = getTotalNumInputChannels();
    const int numOut = getTotalNumOutputChannels();

    for (int ch = numIn; ch < numOut; ++ch)
        buffer.clear(ch, 0, numSamples);

    if (numSamples == 0 || oversampling == nullptr)
        return;

    updateDspParameters(false);

    // Advance block-level smoothers cleanly
    gain1Smoothed.skip(numSamples);
    gain2Smoothed.skip(numSamples);
    trebleSmoothed.skip(numSamples);
    bassSmoothed.skip(numSamples);
    middleSmoothed.skip(numSamples);
    masterSmoothed.skip(numSamples);
    inputVoltsSmoothed.skip(numSamples);

    jcm2203::PreampParameters dspParams;
    dspParams.gain1  = gain1Smoothed.getCurrentValue();
    dspParams.gain2  = gain2Smoothed.getCurrentValue();
    dspParams.high1  = apvts.getRawParameterValue("high1")->load() > 0.5f;
    dspParams.high2  = apvts.getRawParameterValue("high2")->load() > 0.5f;
    dspParams.treble = trebleSmoothed.getCurrentValue();
    dspParams.bass   = bassSmoothed.getCurrentValue();
    dspParams.middle = middleSmoothed.getCurrentValue();
    dspParams.master = masterSmoothed.getCurrentValue();
    dspParams.inputGainVolts = inputVoltsSmoothed.getCurrentValue();

    const int activeCh = std::min(numOut, static_cast<int>(channels.size()));
    doubleScratchBuffer.setSize(activeCh, numSamples, false, false, true);

    for (int ch = 0; ch < activeCh; ++ch)
    {
        const auto* src = buffer.getReadPointer(ch);
        auto* dst = doubleScratchBuffer.getWritePointer(ch);
        for (int i = 0; i < numSamples; ++i)
            dst[i] = static_cast<double>(src[i]);
    }

    juce::dsp::AudioBlock<double> baseBlock(doubleScratchBuffer);
    auto osBlock = oversampling->processSamplesUp(baseBlock);
    const size_t osSamples = osBlock.getNumSamples();

    for (int ch = 0; ch < activeCh; ++ch)
    {
        auto& solver = channels[static_cast<size_t>(ch)];
        solver.setParameters(dspParams);
        double* chPtr = osBlock.getChannelPointer(static_cast<size_t>(ch));
        for (size_t i = 0; i < osSamples; ++i)
            chPtr[i] = solver.processSample(chPtr[i]);
    }

    oversampling->processSamplesDown(baseBlock);

    for (int ch = 0; ch < activeCh; ++ch)
    {
        const auto* src = doubleScratchBuffer.getReadPointer(ch);
        auto* dst = buffer.getWritePointer(ch);
        for (int i = 0; i < numSamples; ++i)
            dst[i] = static_cast<SampleType>(src[i]);
    }
}

juce::AudioProcessorEditor* JCM2203ModdedPreampAudioProcessor::createEditor()
{
    return new juce::GenericAudioProcessorEditor(*this);
}

void JCM2203ModdedPreampAudioProcessor::getStateInformation(juce::MemoryBlock& destData)
{
    if (auto state = apvts.copyState().createXml())
        copyXmlToBinary(*state, destData);
}

void JCM2203ModdedPreampAudioProcessor::setStateInformation(const void* data, int sizeInBytes)
{
    if (auto xml = getXmlFromBinary(data, sizeInBytes))
        if (xml->hasTagName(apvts.state.getType()))
            apvts.replaceState(juce::ValueTree::fromXml(*xml));
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new JCM2203ModdedPreampAudioProcessor();
}
