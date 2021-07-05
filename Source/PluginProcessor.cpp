

#include "PluginProcessor.h"
#include "PluginEditor.h"

#include <Windows.h>
#include <iostream>
#include <sstream>
#include <cmath>

#define DBOUT( s )            \
{                             \
   std::wostringstream os_;    \
   os_ << s;                   \
   OutputDebugStringW( os_.str().c_str() );  \
}

AudioProcessorValueTreeState::ParameterLayout AnotherDelayAudioProcessor::createParameterLayout()
{
	std::vector<std::unique_ptr<RangedAudioParameter>> params;

	{
		using FloatParamPair = std::pair<Identifier, AudioParameterFloat*&>;

		for (auto p : { FloatParamPair(Parameters::gain, gain),
						FloatParamPair(Parameters::delaytime, delaytime),
						FloatParamPair(Parameters::feedback, feedback),
						FloatParamPair(Parameters::mix, mix),
						FloatParamPair(Parameters::lowpass, lowpass),
						FloatParamPair(Parameters::highpass, highpass),
						FloatParamPair(Parameters::flutterfreq, flutterfreq),
						FloatParamPair(Parameters::flutterdepth, flutterdepth),
						FloatParamPair(Parameters::wowfreq, wowfreq),
						FloatParamPair(Parameters::wowdepth, wowdepth),
						FloatParamPair(Parameters::roomsize, roomsize),
						FloatParamPair(Parameters::damping, damping),
						FloatParamPair(Parameters::width, width),
			
			})
		{
			auto& info = Parameters::parameterInfoMap[p.first];
			auto param = std::make_unique<AudioParameterFloat>(p.first.toString(), info.labelName, NormalisableRange<float>(info.lowerLimit, info.upperLimit, info.interval), info.defaultValue);

			p.second = param.get();
			params.push_back(std::move(param));
		}
	}

	{
		auto& info = Parameters::parameterInfoMap[Parameters::reverbenabled];
		auto param = std::make_unique<AudioParameterFloat>(Parameters::reverbenabled.toString(), info.labelName, NormalisableRange<float>(info.lowerLimit, info.upperLimit, info.interval), info.defaultValue);

		reverbenabled = param.get();
		params.push_back(std::move(param));
	}
	return { params.begin(), params.end() };
}


AnotherDelayAudioProcessor::AnotherDelayAudioProcessor()
	: state(*this, nullptr, "PARAMETERS", createParameterLayout()),

#ifndef JucePlugin_PreferredChannelConfigurations
	 AudioProcessor(BusesProperties()
		#if ! JucePlugin_IsMidiEffect
		#if ! JucePlugin_IsSynth
			.withInput("Input", AudioChannelSet::stereo(), true)
		#endif
			.withOutput("Output", AudioChannelSet::stereo(), true)
		#endif
		)
#endif

{
	addParameterListeners();
}

AnotherDelayAudioProcessor::~AnotherDelayAudioProcessor()
{

}


const String AnotherDelayAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool AnotherDelayAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool AnotherDelayAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool AnotherDelayAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double AnotherDelayAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int AnotherDelayAudioProcessor::getNumPrograms()
{
    return 1;   
                
}

int AnotherDelayAudioProcessor::getCurrentProgram()
{
    return 0;
}

void AnotherDelayAudioProcessor::setCurrentProgram (int index)
{
}

const String AnotherDelayAudioProcessor::getProgramName (int index)
{
    return {};
}

void AnotherDelayAudioProcessor::changeProgramName (int index, const String& newName)
{
}


void AnotherDelayAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    
	const int numInputChannels = getTotalNumInputChannels();
	const int delayBufferSize = sampleRate * 10;
	mSampleRate = sampleRate;

	delayBuffer.setSize(numInputChannels, delayBufferSize, false, true);
	wetBuffer.setSize(numInputChannels, samplesPerBlock, false, true);

	dsp::ProcessSpec spec;
	spec.sampleRate = sampleRate;
	spec.maximumBlockSize = samplesPerBlock;
	spec.numChannels = getTotalNumOutputChannels();

	lowPassFilter0.setCoefficients(IIRCoefficients::makeLowPass(44100.0f, 15000.0f));
	lowPassFilter1.setCoefficients(IIRCoefficients::makeLowPass(44100.0f, 15000.0f));
	hiPassFilter0.setCoefficients(IIRCoefficients::makeHighPass(44100.0f, 300.0f));
	hiPassFilter1.setCoefficients(IIRCoefficients::makeHighPass(44100.0f, 300.0f));

	oscFlutterL.setFrequency(1.0);
	oscFlutterL.setSampleRate(44100.0);
	oscFlutterR.setFrequency(1.0);
	oscFlutterR.setSampleRate(44100.0);
	oscWowL.setFrequency(1.0);
	oscWowL.setSampleRate(44100.0);
	oscWowR.setFrequency(1.0);
	oscWowR.setSampleRate(44100.0);

	reverbL.setSampleRate(sampleRate);
	reverbR.setSampleRate(sampleRate);

	updateProcessing();
}

void AnotherDelayAudioProcessor::releaseResources()
{
    
}

void AnotherDelayAudioProcessor::addParameterListeners()
{
	auto& state = getValueTreeState();

	state.addParameterListener(Parameters::gain.toString(), this);
	state.addParameterListener(Parameters::delaytime.toString(), this);
	state.addParameterListener(Parameters::feedback.toString(), this);
	state.addParameterListener(Parameters::mix.toString(), this);
	state.addParameterListener(Parameters::lowpass.toString(), this);
	state.addParameterListener(Parameters::highpass.toString(), this);
	state.addParameterListener(Parameters::flutterfreq.toString(), this);
	state.addParameterListener(Parameters::flutterdepth.toString(), this);
	state.addParameterListener(Parameters::wowfreq.toString(), this);
	state.addParameterListener(Parameters::wowdepth.toString(), this);
	state.addParameterListener(Parameters::reverbenabled.toString(), this);
	state.addParameterListener(Parameters::roomsize.toString(), this);
	state.addParameterListener(Parameters::damping.toString(), this);
	state.addParameterListener(Parameters::width.toString(), this);
}


void AnotherDelayAudioProcessor::parameterChanged(const String& parameterID, float newValue)
{

	if (parameterID == Parameters::gain.toString())
		*gain = newValue;
		updateProcessing();
	if (parameterID == Parameters::delaytime.toString())
		*delaytime = newValue;
		updateProcessing();
	if (parameterID == Parameters::feedback.toString())
		*feedback = newValue;
		updateProcessing();
	if (parameterID == Parameters::mix.toString())
		*mix = newValue;
		updateProcessing();
	if (parameterID == Parameters::lowpass.toString())
		*lowpass = newValue;
		updateProcessing();
	if (parameterID == Parameters::highpass.toString())
		*highpass = newValue;
		updateProcessing();
	if (parameterID == Parameters::flutterfreq.toString())
		*flutterfreq = newValue;
		updateProcessing();
	if (parameterID == Parameters::flutterdepth.toString())
		*flutterdepth = newValue;
		updateProcessing();
	if (parameterID == Parameters::wowfreq.toString())
		*wowfreq = newValue;
		updateProcessing();
	if (parameterID == Parameters::wowdepth.toString())
		*wowdepth = newValue;
		updateProcessing();

	if (parameterID == Parameters::roomsize.toString())
		*roomsize = newValue;
		updateProcessing();
	if (parameterID == Parameters::damping.toString())
		*damping = newValue;
		updateProcessing();
	if (parameterID == Parameters::width.toString())
		*width = newValue;
		updateProcessing();
	if (parameterID == Parameters::reverbenabled.toString())
		*reverbenabled = newValue;
		updateProcessing();
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool AnotherDelayAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    ignoreUnused (layouts);
    return true;
  #else
    
    if (layouts.getMainOutputChannelSet() != AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != AudioChannelSet::stereo())
        return false;

    
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

void AnotherDelayAudioProcessor::updateProcessing()
{
	theDelayEngine.gainInput = (float)(Decibels::decibelsToGain((float)*gain));
	theDelayEngine.delayTimeInput = *delaytime;
	theDelayEngine.feedbackInput = (float)(Decibels::decibelsToGain((float)*feedback));
	theDelayEngine.mixInput = *mix;

	updateOscillator(0);
	updateOscillator(1);

	updateFilter();

	reverbLParameters.roomSize = *roomsize;
	reverbLParameters.damping = *damping;
	reverbLParameters.width = *width;
	reverbL.setParameters(reverbLParameters);

	reverbRParameters.roomSize = *roomsize;
	reverbRParameters.damping = *damping;
	reverbRParameters.width = *width;
	reverbR.setParameters(reverbRParameters);

}

void AnotherDelayAudioProcessor::processBlock(AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
	ScopedNoDenormals noDenormals;
	auto totalNumInputChannels = getTotalNumInputChannels();
	auto totalNumOutputChannels = getTotalNumOutputChannels();

	for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
		buffer.clear(i, 0, buffer.getNumSamples());

	const int delayBufferLength = delayBuffer.getNumSamples();
	const int bufferLength = buffer.getNumSamples();

	for (int channel = 0; channel < totalNumInputChannels; ++channel)
	{
		

		const float* bufferReadPtr = buffer.getReadPointer(channel);
		const float* delayBufferReadPtr = delayBuffer.getReadPointer(channel);
		const float* wetBufferReadPtr = wetBuffer.getReadPointer(channel);
		float* bufferWritePtr = buffer.getWritePointer(channel);
		float* delayBufferWritePtr = delayBuffer.getWritePointer(channel);
		float* wetBufferWritePtr = wetBuffer.getWritePointer(channel);

		
		if (channel == 0)
			fillBuffer(0, bufferLength, delayBufferLength, bufferReadPtr, dBWritePositionL, theDelayEngine.gainInput, lastInputGain);
		else if (channel == 1)
			fillBuffer(1, bufferLength, delayBufferLength, bufferReadPtr, dBWritePositionR, theDelayEngine.gainInput, lastInputGain);
		lastInputGain = theDelayEngine.gainInput;

		float writeValue = {};
		updateFilter(); 					

		
		for (int i = 0; i < bufferLength; i++)
		{
			

			int k;
			float delayTimeInSamples;
			double bpm;

			AudioPlayHead* const ph = getPlayHead();
			AudioPlayHead::CurrentPositionInfo result;
			if (ph != nullptr && ph->getCurrentPosition(result))
				bpm = result.bpm;

			
			if (channel == 0)
				k = dBWritePositionL;
			else if (channel == 1)
				k = dBWritePositionR;

			
			float delayTimeInput = (120000.0 / theDelayEngine.delayTimeInput) / bpm; 

			delayTimeInSamples = ((float)mSampleRate * delayTimeInput / 1000.0);
			float delaySampleFloor = (floor)(delayTimeInSamples);

			
			float cleansig = *bufferWritePtr;

			if (delayTimeInSamples == 0.f)
				*wetBufferWritePtr = cleansig;
			else 
			{
				updateOscillator(channel); 

				if (delayTimeInSamples - delaySampleFloor != 0) {
					delayTimeInSamples -= (delayTimeInSamples - delaySampleFloor);
				}
				
				writeValue = calcWriteValue(channel, delayBuffer, k, i, delayBufferLength, delayTimeInSamples, mod);
				*wetBufferWritePtr = writeValue * sqrt(theDelayEngine.mixInput);
				
				if (channel == 0)
				{
					*wetBufferWritePtr = lowPassFilter0.processSingleSampleRaw(*wetBufferWritePtr);
					*wetBufferWritePtr = hiPassFilter0.processSingleSampleRaw(*wetBufferWritePtr);
				}
				else if (channel == 1) {
					*wetBufferWritePtr = lowPassFilter1.processSingleSampleRaw(*wetBufferWritePtr);
					*wetBufferWritePtr = hiPassFilter1.processSingleSampleRaw(*wetBufferWritePtr);		
				}
				
				*wetBufferWritePtr = (1 / atan(2)) * atan(2 * *wetBufferWritePtr);

				
				if (*reverbenabled == 1)
				{
					if (channel == 0)
						reverbL.processMono(wetBufferWritePtr, 1);
					else if (channel == 1)
						reverbR.processMono(wetBufferWritePtr, 1);
				}
				
				*bufferWritePtr = *wetBufferWritePtr + (sqrt(1 - theDelayEngine.mixInput) * cleansig);

				
				delayBuffer.addFromWithRamp(channel, (k + i) % delayBufferLength, bufferWritePtr, 1, theDelayEngine.feedbackInput, lastFeedbackGain);
				lastFeedbackGain = theDelayEngine.feedbackInput;
			}
			bufferWritePtr++;
			wetBufferWritePtr++;	
		
		}
	

		if (channel == 0)
		{
			dBWritePositionL += bufferLength;
			dBWritePositionL %= delayBufferLength;
		}
		else if (channel == 1) {
			dBWritePositionR += bufferLength;
			dBWritePositionR %= delayBufferLength;
		}
	}
}



float AnotherDelayAudioProcessor::calcWriteValue(int channel, AudioBuffer<float>& buffer, int k, int i, int delayBufferLength, float delayTimeInSamples, float mod)
{
	float kk = (float)k;
	float ii = (float)i;
	float delayBufferLengthF = (float)delayBufferLength;

	float floorValue0 = buffer.getSample(channel, ((int)(floor)(kk + ii - delayTimeInSamples - mod + delayBufferLengthF) - 1) % delayBufferLength);
	float floorValue = buffer.getSample(channel, ((int)(floor)(kk + ii - delayTimeInSamples - mod + delayBufferLengthF) - 0) % delayBufferLength);
	float floorValue1 = buffer.getSample(channel, ((int)(ceil)(kk + ii - delayTimeInSamples - mod + delayBufferLengthF) + 0) % delayBufferLength);
	float floorValue2 = buffer.getSample(channel, ((int)(ceil)(kk + ii - delayTimeInSamples - mod + delayBufferLengthF) + 1) % delayBufferLength);


	return interpolate(floorValue0, floorValue, floorValue1, floorValue2, delayTimeInSamples, mod);
}


void AnotherDelayAudioProcessor::fillBuffer(int channel, const int bufferLength, const int delayBufferLength, const float* bufferReadPtr,
	int dBWritePosition, float startGain, float endGain)
{
	if (delayBufferLength > bufferLength + dBWritePosition)
	{
		delayBuffer.copyFromWithRamp(channel, dBWritePosition, bufferReadPtr, bufferLength, startGain, endGain);
	}
	else {
		const int bufferRemaining = delayBufferLength - dBWritePosition;
		const float midGain = lastInputGain + ((endGain - startGain) / bufferLength) * (bufferRemaining / bufferLength);
		delayBuffer.copyFromWithRamp(channel, dBWritePosition, bufferReadPtr, bufferRemaining, startGain, midGain);
		delayBuffer.copyFromWithRamp(channel, 0, bufferReadPtr + bufferRemaining, bufferLength - bufferRemaining, midGain, endGain);
	}
}


void AnotherDelayAudioProcessor::fetchDelay(AudioBuffer<float>& buffer, int channel, const int feedbackBufferLength,
	const int delayBufferLength, const float* feedbackBufferPtr, const float* delayBufferPtr, float startGain, float endGain)
{
	int delayTimeInput = *delaytime;
	int delayTimeInSamples = (mSampleRate * delayTimeInput / 1000.0);
	const int readPosition = (int)(delayBufferLength + dBWritePositionL - delayTimeInSamples) % delayBufferLength;

	if (delayBufferLength > feedbackBufferLength + readPosition)
	{
		buffer.copyFrom(channel, 0, delayBufferPtr + readPosition, feedbackBufferLength);
	}
	else {
		const int bufferRemaining = delayBufferLength - readPosition;
		buffer.copyFrom(channel, 0, delayBufferPtr + readPosition, bufferRemaining);
		buffer.copyFrom(channel, bufferRemaining, delayBufferPtr, feedbackBufferLength - bufferRemaining);
	}
}


void AnotherDelayAudioProcessor::sendFeedback(AudioBuffer<float>& buffer, int channel, const int feedbackBufferLength, const int delayBufferLength, float* feedbackBufferWritePtr,
	float startGain, float endGain)
{
	if (delayBufferLength > feedbackBufferLength + dBWritePositionL)
	{
		delayBuffer.addFromWithRamp(channel, dBWritePositionL, feedbackBufferWritePtr, feedbackBufferLength, startGain, endGain);
	}
	else {
		const float bufferRemaining = delayBufferLength - dBWritePositionL;
		const float midGain = lastInputGain + ((endGain - startGain) / feedbackBufferLength) * (bufferRemaining / feedbackBufferLength);
		delayBuffer.addFromWithRamp(channel, dBWritePositionL, feedbackBufferWritePtr, bufferRemaining, startGain, midGain);
		delayBuffer.addFromWithRamp(channel, 0, feedbackBufferWritePtr, feedbackBufferLength - bufferRemaining, midGain, endGain);
	}
}


AudioProcessorValueTreeState& AnotherDelayAudioProcessor::getValueTreeState()
{
	return state;
}

bool AnotherDelayAudioProcessor::hasEditor() const
{
	return true; 
}

AudioProcessorEditor* AnotherDelayAudioProcessor::createEditor()
{
    return new AnotherDelayAudioProcessorEditor (*this);
}

void AnotherDelayAudioProcessor::updateFilter()
{
	

	lowPassFilter0.setCoefficients(IIRCoefficients::makeLowPass(44100.0f, *lowpass));
	lowPassFilter1.setCoefficients(IIRCoefficients::makeLowPass(44100.0f, *lowpass));
	hiPassFilter0.setCoefficients(IIRCoefficients::makeHighPass(44100.0f, *highpass));
	hiPassFilter1.setCoefficients(IIRCoefficients::makeHighPass(44100.0f, *highpass));
}

void AnotherDelayAudioProcessor::getStateInformation (MemoryBlock& destData)
{
	MemoryOutputStream stream(destData, false);
	state.state.writeToStream(stream);
}

void AnotherDelayAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
	ValueTree tree = ValueTree::readFromData(data, sizeInBytes);
	if (tree.isValid()) {
		state.state = tree;
	}
}

double AnotherDelayAudioProcessor::updateOscillator(int channel)
{
	auto target = 0.0;
	auto flutterDepth = *flutterdepth * 0.5;
	auto wowDepth = *wowdepth * 0.5;
	auto totalDepth = flutterDepth + wowDepth;
	auto followSpeed = 1.0;

	if (channel == 0)
	{
		auto oscFlutterLValue = oscFlutterL.nextSample(*flutterfreq);
		if (flutterDepth < 0.0)
		{
			flutterDepth *= -1.0;
			oscFlutterLValue = 1.0 - oscFlutterLValue;
		}
		target += flutterDepth * oscFlutterLValue;

		auto oscWowLValue = oscWowL.nextSample(*wowfreq);
		if (wowDepth < 0.0)
		{
			wowDepth *= -1.0;
			oscWowLValue = 1.0 - oscWowLValue;
		}
		target += wowDepth * oscWowLValue;
		
		if (flutterDepth && wowDepth == 0.0)
			auto followSpeed = 10.0;
		else
			auto followSpeed = 1.0;
	
		oscLPosition += (target - oscLPosition) * followSpeed * (double)(1.0 / mSampleRate);
		mod = oscLPosition * (double)mSampleRate;
	}
	else if (channel == 1)
	{ 
		auto oscFlutterRValue = oscFlutterR.nextSample(*flutterfreq);
		if (flutterDepth < 0.0)
		{
			flutterDepth *= -1.0;
			oscFlutterRValue = 1.0 - oscFlutterRValue;
		}
		target += flutterDepth * oscFlutterRValue;

		auto oscWowRValue = oscWowR.nextSample(*wowfreq);
		if (wowDepth < 0.0)
		{
			wowDepth *= -1.0;
			oscWowRValue = 1.0 - oscWowRValue;
		}
		target += wowDepth * oscWowRValue;

		if (flutterDepth && wowDepth == 0.0)
			auto followSpeed = 10.0;
		else
			auto followSpeed = 1.0;

		oscRPosition += (target - oscRPosition) * followSpeed * (double)(1.0 / mSampleRate);
		mod = (oscRPosition * (double)mSampleRate);
	}
	return mod;


	
}




AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new AnotherDelayAudioProcessor();
}


