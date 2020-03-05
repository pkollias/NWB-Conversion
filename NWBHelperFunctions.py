from os import getcwd, path
import numpy as np
from scipy.signal import savgol_filter

import pynwb
from pynwb import NWBHDF5IO, ProcessingModule
from pynwb.base import Images
from pynwb.image import ImageSeries

from matplotlib import pyplot as plt

def OpenNWBFile(nwb_file_name, mode='a', verbose=False):
    """
        Opens a Neurodata Without Borders (NWB) file.

        Usage: nwb_file, nwb_io = OpenNWBFile(nwb_file_name)

        :param nwb_file_name: {str} Name of NWB file to open.
        :param mode: [optional] {str} Mode of opening. Default is appending (mode='a'). Can be 'r' or 'w' for reading and writing (will overwrite)
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :return: {pynwb.file.NWBFile} Handle for NWB data file
        :return: {pynwb.NWBHDF5IO} Handle for NWB data file IO
            """

    # Check to see if user passed a filename or file
    if not nwb_file_name:
        raise ValueError("Must pass a NWB file name to open.")
    # Check to make sure file exists
    if not path.isfile(nwb_file_name):
        raise FileNotFoundError("NWB file doesn't exist.")
    # Initialize NWB file
    if verbose: print("Opening NWB file %s..." % (nwb_file_name), end='')
    nwb_io = NWBHDF5IO(nwb_file_name, mode=mode)
    nwb_file = nwb_io.read()
    if verbose: print("done.")

    return nwb_file, nwb_io

def CloseNWBFile(nwb_file, nwb_io, verbose=False):
    """
        Closes a Neurodata Without Borders (NWB) file.

        Usage: CloseNWBFile(nwb_file, nwb_io)

        :param nwb_file: {NWBFile} NWB file.
        :param nwb_io: {NWBHDF5IO} File handle for writing NWB file.
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
            """

    # Write to file and close
    if verbose: print("Writing NWB file and closing.")
    nwb_io.write(nwb_file)
    nwb_io.close()


def AddStimulusImageToNWB(nwb_file, image_info, verbose=False):
    """
        Adds stimulus images to Neurodata Without Borders (NWB) file.

        :param nwb_file: {NWB file or str} Either pointer to existing NWB file or a string so NWB file can be opened.
        :param image_info: {dict} Dictionary of images to import.
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :return: {str} filename of NWB file (empty if error)
        """

    # Check to see if user passed a filename or file
    if not nwb_file:
        raise ValueError("Must pass a NWB file.")
    if type(nwb_file) == pynwb.file.NWBFile:
        if verbose: print("Using passed NWB file.")
    elif type(nwb_file) == str:
        nwb_file_name = nwb_file
        [nwb_file, nwb_io] = OpenNWBFile(nwb_file_name, verbose=verbose)
    else:
        raise ValueError("NWB file was not valid.")

    # Check to see if user passed image files
    if not image_info:
        if nwb_io: nwb_io.close()
        raise ValueError("Must pass image information to add to NWB file.")

    # Make sure behavior module exists in the file
    if 'behavior' not in nwb_file.processing.keys():
        if verbose: print("Behavior processing module does not exist. Creating.")
        signal_module = ProcessingModule(name='behavior',
                                         description="Processing module for behavior signals.")
        nwb_file.add_processing_module(signal_module)

    # Create interface for the stimulus images
    if verbose: print("Creating interface for holding stimulus images.")
    stim_img_interface = Images("Stimulus Images", description="Interface for storing the raw stimulus images used in experiment.")

    # Loop through and add each image
    for cur_image_ind, cur_image_info in enumerate(image_info):
        # Load image from file
        if not path.isfile(cur_image_info['file']):
            raise FileNotFoundError("Image file '%s' doesn't exist." % (cur_image_info['file']))
        img = plt.imread(cur_image_info['file'])

        # Save image to interface
        stim_img_interface.create_image(name=cur_image_info['name'],
                                        data=img,
                                        resolution=cur_image_info['resolution'],
                                        description=cur_image_info['description'])
        if verbose: print("\tAdded image %s to interface." % (cur_image_info['name']))

    # Add image interface to behavior module in NWB file
    if verbose: print("Adding stimulus image data interface to behavior module.")
    nwb_file.processing['behavior'].add(stim_img_interface)

    # Write the file
    if nwb_io:
        if verbose: print("Writing NWB file and closing.")
        nwb_io.write(nwb_file)
        nwb_io.close()
        return nwb_file_name
    else:
        if verbose: print("Returning NWB file variable")
        return nwb_file


def AddTrialDataToNWB(nwb_file, trial_data, trial_descriptions=None, verbose=False):
    """
        Adds information about a trial (or trials) to Neurodata Without Borders (NWB) file.

        :param nwb_file: {NWB file or str} Either pointer to existing NWB file or a string so NWB file can be opened.
        :param trial_data: {list} List of trial information (stored in a dict) to be put in NWB file.
        :param trial_descriptions: {dict} A dictionary carrying descriptions of the columns for trials. Helps make NWB file more user-readable.
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :return: {str} filename of NWB file (empty if error)
        """

    # Check to see if user passed a filename or file
    if not nwb_file:
        raise ValueError("Must pass a NWB file.")
    if type(nwb_file) == pynwb.file.NWBFile:
        if verbose: print("Using passed NWB file.")
    elif type(nwb_file) == str:
        nwb_file_name = nwb_file
        [nwb_file, nwb_io] = OpenNWBFile(nwb_file_name, verbose=verbose)
    else:
        raise ValueError("NWB file was not valid.")

    # Verify the trial information passed
    if not trial_data:
        raise ValueError("Must pass trial data to add to NWB file.")
    if type(trial_data) is not list:
        raise ValueError("Trial data should be in the format of a list of dictionaries.")

    # Process keys in input
    for cur_key in trial_data[0].keys():
        # Don't process the start and stop times
        if (cur_key == 'start_time') or (cur_key == 'stop_time'):
            continue

        # Create NWB columns if necessary
        if (nwb_file.trials is None) or (cur_key not in nwb_file.trials.colnames):
            #Check to see if description of column was provided
            if cur_key in trial_descriptions.keys():
                cur_key_desc = trial_descriptions[cur_key]
            else:
                cur_key_desc = 'Not provided.'
            # Add column to NWB file
            if verbose: print("Adding column %s to trial table. Description: %s" % (cur_key, cur_key_desc))
            nwb_file.add_trial_column(name=cur_key, description=cur_key_desc)

        # Process inputs to make sure they are all the same length and also check if we have to handle strings
        max_len = 0
        data_is_str = False
        f_vect_is_str = np.vectorize(lambda x: isinstance(x, str))
        for cur_trial in trial_data:
            if hasattr(cur_trial[cur_key], '__iter__'):
                max_len = np.max((max_len, len(cur_trial[cur_key])))
            else:
                max_len = np.max((max_len, 1))
            data_is_str = data_is_str or np.any(f_vect_is_str(np.array(cur_trial[cur_key])))

        for cur_trial in trial_data:
            if data_is_str:
                # We have to serialize strings for NWB to handle them.  Seems to have issues with arrays of strings.
                temp_data = ''
                for cur_val_ind, cur_val in enumerate(cur_trial[cur_key]):
                    temp_data = temp_data + str(cur_val) + ';'
                cur_trial[cur_key] = temp_data
            elif (max_len > 1):
                temp_data = np.full((max_len,), -99, dtype=type(cur_trial[cur_key][0]))
                temp_data[0:len(cur_trial[cur_key])] = cur_trial[cur_key]
                cur_trial[cur_key] = temp_data


    # Loop through trials, adding each one
    if verbose: print("Adding trial data to NWB file.")
    for cur_trial_ind, cur_trial in enumerate(trial_data):
        nwb_file.add_trial(**cur_trial)

    # Write the file
    if nwb_io:
        if verbose: print("Writing NWB file and closing.")
        nwb_io.write(nwb_file)
        nwb_io.close()
        return nwb_file_name
    else:
        if verbose: print("Returning NWB file variable")
        return nwb_file


def AddMesoscaleImagingDataToNWB(nwb_file, tiff_files=None, img_info=None, timestamps=None, load_data=False, verbose=False):
    """
        Adds information about a trial (or trials) to Neurodata Without Borders (NWB) file.

        :param nwb_file: {NWB file or str} Either pointer to existing NWB file or a string so NWB file can be opened.
        :param tiff_files: {list} List of TIFF files. All will be loaded (in order provided)
        :param img_info: {dict} A dictionary carrying information about the images.
        :param load_data: [optional] {bool} whether to load images from TIFF stack or leave it as a reference. Default is false.
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :return: {str} filename of NWB file (empty if error)
        """

    # Check to see if user passed a filename or file
    nwb_io = None
    if not nwb_file:
        raise ValueError("Must pass a NWB file.")
    if type(nwb_file) == pynwb.file.NWBFile:
        if verbose: print("Using passed NWB file.")
    elif type(nwb_file) == str:
        nwb_file_name = nwb_file
        [nwb_file, nwb_io] = OpenNWBFile(nwb_file_name, verbose=verbose)
    else:
        raise ValueError("NWB file was not valid.")

    # Define default image information
    if verbose: print("Parsing image information.")
    def_image_info = {'resolution': -1.0,
                      'bits_per_pixel': 0,
                      'starting_time': -1.0,
                      'rate': -1.0,
                      'comments': 'No comments provided.',
                      'description': 'No description provided.'}
    # Update the user-specified keys
    if img_info is not None:
        def_image_info.update(img_info)
    img_info = def_image_info

    # Check whether we were passed timestamps for the images
    if (timestamps is None) or (len(timestamps) == 0):
        # Use starting time and rate
        starting_time = img_info['starting_time']
        rate = img_info['rate']
        timestamps = None
    else:
        # Use timestamps
        starting_time = None
        rate = None

    # Load data?
    if load_data:
        if verbose: print("Loading data from files to store in NWB.")
        # Load all of the data from file
        img = []
        for cur_img_file in tiff_files:
            # Read with the PIL library because something is wrong with header for TIFFfile
            cur_img = io.imread(cur_img_file, plugin='pil')
            # Append to overall image
            if len(img) == 0:
                img = cur_img
            else:
                img = np.concatenate([img, cur_img], axis=0)
            if verbose: print("\tLoaded %s (contains %d images of %d x %d)." % (cur_img_file, cur_img.shape[0], cur_img.shape[1], cur_img.shape[2]))

        if verbose: print("Creating image series.")
        image_series = ImageSeries("Raw Mesoscale Images from TIFF Stack",
                                   data=img,
                                   unit='fluorescence',
                                   format='TIFF',
                                   bits_per_pixel=img_info['bits_per_pixel'],
                                   dimension=[img.shape[1], img.shape[2]],
                                   resolution=img_info['resolution'],
                                   conversion=1.0,
                                   timestamps=timestamps,
                                   starting_time=starting_time,
                                   rate=rate,
                                   comments=img_info['comments'] + " Images from files: %s" % ([f for f in tiff_files]),
                                   description=img_info['description'])
    else: # don't load data, store by reference
        # Need to grab the size of each image stack so we can keep a store of starting image frames
        num_frames = []
        for cur_img_file, cur_img_file_ind in enumerate(tiff_files):
            # Read with the PIL library because something is wrong with header for TIFFfile
            cur_img.clear()
            cur_img = io.imread(cur_img_file, plugin='pil')
            num_frames.append(cur_img.shape[0])
        # Pop off the last one and push on a first one
        num_frames.insert(0, 0)
        num_frames.pop()
        starting_frames = np.array(num_frames).cumsum(axis=0)

        image_series = ImageSeries("Raw Mesoscale Images from TIFF Stack",
                                   data=img,
                                   unit='fluorescence',
                                   format='external',
                                   external_file=tiff_files,
                                   starting_frame=starting_frames,
                                   bits_per_pixel=img_info['bits_per_pixel'],
                                   dimension=[cur_img.shape[1], cur_img.shape[2]],
                                   resolution=img_info['resolution'],
                                   conversion=1.0,
                                   timestamps=timestamps,
                                   starting_time=starting_time,
                                   rate=rate,
                                   comments=img_info['comments'] + " Images from files: %s" % ([f for f in tiff_files]),
                                   description=img_info['description'])

    # Add to the file
    nwb_file.add_acquisition(image_series)

    # Write the file
    if nwb_io:
        if verbose: print("Writing NWB file and closing.")
        nwb_io.write(nwb_file)
        nwb_io.close()
        return nwb_file_name
    else:
        if verbose: print("Returning NWB file variable")
        return nwb_file

def ExtractTimestampsFromAnalogSignal(analog_signal, threshold, smooth=None, index=None, verbose=False):
    """
        Extracts threshold crossings from an analog signal. Can be useful for extracting event times.

        :param analog_signal: {pynwb.base.TimeSeries} A TimeSeries (or subclass) pulled from a NWB file.
        :param threshold: {float or int} Threshold to use to mark events in the analog signal.
        :param smooth: [optional] {dict} If specified, smooths analog signal according to dictionary. Default is None.
        :param index: [optional] {int} If the analog signal has multiple channels, you must specify which to use.
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :return: {array, array} onset and offset times (in seconds) of events
        """

    # Check to see if user passed a valid timeseries
    if (not analog_signal) or not issubclass(type(analog_signal), pynwb.base.TimeSeries):
        raise ValueError("Must pass an TimeSeries as the analog signal")
    # Check to see if threshold is valid
    if not threshold:
        raise ValueError("Must pass a threshold.")

    # Extract the data for the current signal index
    if not index:
        if analog_signal.data.shape[1] > 1:
            raise ValueError("Can only threshold one signal at a time. If passing a analog signal with multiple channels, must specify the index.")
        index = 0
    cur_data = analog_signal.data[:, index]

    # Create a time series for the analog signal
    if not analog_signal.timestamps:
        if verbose: print("Creating timestamps for analog signal (starting at %f seconds with rate %f)." % (analog_signal.starting_time, analog_signal.rate))
        cur_ts = analog_signal.starting_time + np.arange(0, (cur_data.shape[0]/analog_signal.rate), (1/analog_signal.rate))
    else:
        if verbose: print("Taking timestamps from analog signal.")
        cur_ts = analog_signal.timestamps

    # Smooth the signal
    if smooth:
        if smooth['type'].lower() == 'boxcar':
            #Create boxcar and smooth
            if verbose: print("Smoothing with a boxcar of width %d." % (smooth['width']))
            box_car = np.ones(smooth['width'])/smooth['width']
            cur_data = np.convolve(cur_data, box_car, mode='valid')
        elif smooth['type'].lower() in ['hanning', 'hamming', 'bartlett', 'blackman']:
            # Use built-in numpy functions to do other smoothing.
            if verbose: print("Smoothing with a %s kernel of width %d." % (smooth['type'], smooth['width']))
            kernel = eval("np.%s(smooth['width'])" % (smooth['type'].lower()))
            cur_data = np.convolve(cur_data, kernel, mode='valid')
        elif smooth['type'].lower() == 'savgol_filter':
            cur_data = savgol_filter(cur_data, smooth['width'], smooth['order'])

    # Threshold the signal
    if verbose: print("Thresholding signal at %f." % (threshold))
    thresh_signal = (np.array(cur_data) >= threshold).astype(int)

    # Differentiate to get transitions to up -- pad with a 0 at the beginning and the end
    if verbose: print("Differentiating signal to find transitions")
    diff_signal = np.ediff1d(np.concatenate([np.zeros((1,)), thresh_signal, np.zeros((1,))], axis=0))

    # Onset times
    onset_index = np.where(diff_signal == 1)[0]
    onset_times = cur_ts[onset_index]
    if verbose: print("Found %d onsets." % (len(onset_index)))

    # Offset times
    offset_index = np.where(diff_signal == -1)[0]
    offset_times = cur_ts[offset_index]
    if verbose: print("Found %d offsets." % (len(offset_index)))

    # Plot example trace with threshold
    if verbose:
        # Create figure, plot the start and end of the signal with threshold and markers
        fig = plt.figure()
        axes = fig.subplots(1, 2)
        fig.suptitle('Analog Signal %s w/ Threshold' % (analog_signal.name))

        #Plot the start of the signal
        plt_start = int(np.max([0, round(onset_index[0] - 0.1 * analog_signal.rate)]))
        plt_end = int(np.min([cur_data.shape[0], round(plt_start + 1 * analog_signal.rate)]))
        # Plot the raw data
        axes[0].plot(cur_ts[plt_start:plt_end], cur_data[plt_start:plt_end], label=analog_signal.name)
        # Plot the threshold
        axes[0].plot(axes[0].get_xlim(), threshold*np.ones(2,), label='Threshold')
        # Plot the onset and offset times
        onset_index_good = tuple(np.logical_and([onset_index <= plt_end], [onset_index >= plt_start]))
        axes[0].plot(onset_times[onset_index_good], cur_data[onset_index[onset_index_good]], 'g*', label='Onsets')
        offset_index_good = tuple(np.logical_and([offset_index <= plt_end], [offset_index >= plt_start]))
        axes[0].plot(offset_times[offset_index_good], cur_data[offset_index[offset_index_good]], 'ro', label='Offsets')
        # Labels, legend and show
        axes[0].set_xlabel('Time (s)')
        axes[0].set_ylabel('Signal')
        axes[0].legend()

        # Plot the end of the signal
        plt_end = int(np.min([cur_data.shape[0], round(onset_index[-1] + 0.05 * analog_signal.rate)]))
        plt_start = int(np.max([0, round(plt_end - 0.15 * analog_signal.rate)]))
        # Plot the raw data
        axes[1].plot(cur_ts[plt_start:plt_end], cur_data[plt_start:plt_end], label=analog_signal.name)
        # Plot the threshold
        axes[1].plot(axes[1].get_xlim(), threshold * np.ones(2, ), label='Threshold')
        # Plot the onset and offset times
        onset_index_good = tuple(np.logical_and([onset_index <= plt_end], [onset_index >= plt_start]))
        axes[1].plot(onset_times[onset_index_good], cur_data[onset_index[onset_index_good]], 'g*', label='Onsets')
        offset_index_good = tuple(np.logical_and([offset_index <= plt_end], [offset_index >= plt_start]))
        axes[1].plot(offset_times[offset_index_good], cur_data[offset_index[offset_index_good]], 'ro', label='Offsets')
        # Labels, legend and show
        axes[1].set_xlabel('Time (s)')
        axes[1].set_ylabel('Signal')
        axes[1].legend()

        plt.show()

    return onset_times, offset_times
