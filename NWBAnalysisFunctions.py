from os import path
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

from NWBHelperFunctions import OpenNWBFile
from pynwb import file, NWBFile, NWBHDF5IO, ProcessingModule
from pynwb.misc import DecompositionSeries
from pynwb.ecephys import ElectricalSeries, LFP



def NWBFilterLFP(nwb_file, signal_name=None, elec_ids=None, filter_info=None, downsample=None, copy_file=False, verbose=True):
    """
    Does a wavelet transform on electrophysiology data from NWB file. Saves in a new processing module
    in the same NWB file.

    :param nwb_file: {NWB file or str} Either a handle to the open NWB file or a string with the file name.
    :param signal_name: {str} Name of the signal to process. This should be the name of the ElectricalSeries.  Will search for this in raw or LFP module.
    :param elec_ids: [optional] {List of integers} List of electrodes to decompose. Default is all electrodes.
    :param filter_info: [optional] {dict} Information about filter to apply to data. Should contain 'name' of filter (butter, cheby1, cheby2, ellip), 'filt_order' specifies the filter order, and 'freq_range' specifies the band-pass frequencies. Other filters need passband or stopband dBs specified
    :param downsample: [optional] {int} If specified, downsample signal by this ratio. Default is None (keep original)
    :param copy_file: [optional] {bool} Whether to copy the file (default is False). If true, then a shallow copy of the file is created.
    :param verbose: [optional] {bool} Whether to print updates of progress. Default is False.
    :return: {NWB file or str} Depending on type of nwb_file input passed, either returns an open NWB file handle or the string of the file name.
    """

    # Check to see if user passed a filename or file
    if not nwb_file:
        raise ValueError("Must pass a NWB file.")
    if type(nwb_file) == file.NWBFile:
        if copy_file:
            raise ValueError("If copying file, you must pass file name (need the original IO).")
        if verbose: print("Using passed NWB file.")
    elif type(nwb_file) == str:
        nwb_file_name = nwb_file
        [nwb_file, nwb_io] = OpenNWBFile(nwb_file_name, verbose=verbose)
    else:
        raise ValueError("NWB file was not valid.")

    # Check to see if we should do a copy of the file first
    if copy_file:
        nwb_file = nwb_file.copy()

    # Parse name of signal we want to transform, grab electrical series
    if signal_name in nwb_file.acquisition.keys():
        # Then use raw data
        cur_es = nwb_file.acquisition[signal_name]
        cur_es_data = nwb_file.acquisition[signal_name].data
        es_elec_ids = nwb_file.acquisition[signal_name].electrodes.table[:,0]
        cur_es_rate = nwb_file.acquisition[signal_name].rate
        cur_es_starting_time = nwb_file.acquisition[signal_name].starting_time
    elif ('ecephys' in nwb_file.processing.keys()) and (signal_name in nwb_file.processing['ecephys'].data_interfaces.keys()):
        # Then use LFP data from this interface -- concatenate them all
        cur_es = nwb_file.processing['ecephys'][signal_name].electrical_series
        cur_es_data = None
        for cur_key in cur_es.keys():
            if cur_es_data is None:
                cur_es_data = cur_es[cur_key].data
                es_elec_ids = cur_es[cur_key].electrodes.table[:,0]
                cur_es_rate = cur_es[cur_key].rate
                cur_es_starting_time = cur_es[cur_key].starting_time
            else:
                cur_es_data = np.concatenate((cur_es_data, cur_es[cur_key].data), axis=1)
                es_elec_ids = np.concatenate((es_elec_ids, cur_es[cur_key].electrodes.table[:,0]))
    else:
        raise ValueError("Couldn't find a signal named '%s'." % (signal_name))

    # Create new processing module for transformed data
    # Make sure extracellular ephys module exists in the file
    if 'ecephys' not in nwb_file.processing.keys():
        if verbose: print("Electrophysiology processing module does not exist. Creating.")
        signal_module = ProcessingModule(name='ecephys',
                                         description="Processing module for electrophysiology signals.")
        nwb_file.add_processing_module(signal_module)

    # Check our list of electrode ids
    if elec_ids is None:
        elec_ids = es_elec_ids
    if type(elec_ids) is not np.ndarray:
        elec_ids = np.array(elec_ids)
    elec_ids = elec_ids[np.in1d(elec_ids, es_elec_ids)]

    # Limit our data to the requested electrodes
    cur_es_data = cur_es_data[:, np.in1d(es_elec_ids, elec_ids)]
    cur_es_elec_ids = es_elec_ids[np.in1d(es_elec_ids, elec_ids)]

    # Create electrode table for new electrodes
    elec_table_ind = np.ones((len(cur_es_elec_ids), 1)) * np.NaN
    for cur_elec_ind, cur_elec in enumerate(cur_es_elec_ids):
        # In order to create electrode table, we have to have indexes for each channel into electrode list in NWB file
        cur_elec_table_ind = np.where(nwb_file.electrodes[:, 0] == cur_elec)[0]
        if len(cur_elec_table_ind) == 0:
            raise ValueError("Couldn't find electrode %d in NWB file list." % (cur_elec))
        elec_table_ind[cur_elec_ind] = cur_elec_table_ind
    elec_table_ind = elec_table_ind.transpose().tolist()[0]  # Convert to list for create_electrode_table_region
    cur_es_electrode_table = nwb_file.create_electrode_table_region(elec_table_ind, "Electrodes %s" % (elec_ids))

    # Parse filter type
    if filter_info is None:
        filt_str = "Butterworth: Order %d, %2.1f to %2.1f Hz" % (3, 0.5, 120)
        sos = signal.butter(3, np.array([0.5, 120]) / cur_es_rate, btype='bandpass', output='sos')
    elif filter_info['name'].lower() == 'butter':
        filt_str = "Butterworth: Order %d, %2.1f to %2.1f Hz" % (filter_info['filt_order'], filter_info['freq_range'][0], filter_info['freq_range'][1])
        sos = signal.butter(filter_info['filt_order'], np.array(filter_info['freq_range']) / cur_es_rate, btype='bandpass', output='sos')
    elif filter_info['name'].lower() == 'cheby1':
        filt_str = "Chebyshov Type-I Filter: Order %d, Bandpass Ripple %3.2f dB, %2.1f to %2.1f Hz" % \
                          (filter_info['filt_order'], filter_info['ripple_dB'], filter_info['freq_range'][0], filter_info['freq_range'][1])
        sos = signal.cheby1(filter_info['filt_order'], filter_info['ripple_dB'], np.array([0.5, 120]) / cur_es_rate, btype='bandpass', output='sos')
    elif filter_info['name'].lower() == 'cheby2':
        filt_str = "Chebyshov Type-II Filter: Order %d, Bandpass Ripple %3.2f dB, %2.1f to %2.1f Hz" % \
                          (filter_info['filt_order'], filter_info['ripple_dB'], filter_info['freq_range'][0], filter_info['freq_range'][1])
        sos = signal.cheby2(filter_info['filt_order'], filter_info['ripple_dB'], np.array([0.5, 120]) / cur_es_rate, btype='bandpass', output='sos')
    elif filter_info['name'].lower() == 'ellip':
        filt_str = "Elliptic Filter: Order %d, Bandpass Ripple %3.2f dB, Stopband %3.2f dB, %2.1f to %2.1f Hz" % \
                          (filter_info['filt_order'], filter_info['ripple_dB'], filter_info['stop_dB'], filter_info['freq_range'][0], filter_info['freq_range'][1])
        sos = signal.ellip(filter_info['filt_order'], filter_info['ripple_dB'], filter_info['stop_dB'], np.array([0.5, 120]) / cur_es_rate, btype='bandpass', output='sos')
    if verbose: print("Created filter: %s..." % (filt_str))

    # Apply filter to original data
    if verbose: print("Applying filter to the original data...")
    lfp_data = signal.sosfiltfilt(sos, cur_es_data, axis=0)
    lfp_rate = cur_es.rate

    # Decimate if requested
    if downsample is not None:
        # Downsample the current data (note: this applies an extra filter)
        if verbose: print("Decimating data by a factor of %d..." % (downsample))
        cur_data = signal.decimate(lfp_data, downsample, axis=0, zero_phase=True)
        lfp_rate = cur_es.rate/downsample

    # Create electrical series for filtered data
    lfp_es = ElectricalSeries(name='LFP, filtered from ' + signal_name,
                              data=lfp_data,
                              electrodes=cur_es_electrode_table,
                              starting_time=cur_es.starting_time,
                              rate=lfp_rate,
                              resolution=cur_es.resolution,
                              conversion=cur_es.conversion,
                              # this is to convert the raw data (in mV) to V as expected by NWB
                              comments=cur_es.comments + "; LFPs filtered with: " + filt_str,
                              description="LFP data, filtered from " + signal_name)

    # Create lfp data interface
    cur_lfp = LFP(lfp_es, name='LFP')

    # Add data interface to behavior module in NWB file
    if verbose: print("Adding LFP to electrophysiology module.")
    nwb_file.processing['ecephys'].add(cur_lfp)

    # Write the file
    if nwb_io is not None:
        if copy_file:
            nwb_copy_file_name = path.splitext(nwb_file_name)[0] + '_LFPFiltered' + path.splitext(nwb_file_name)[1]
            if verbose: print("Did a shallow copy of NWB file, so saving to new file: %s" % (nwb_copy_file_name))
            with NWBHDF5IO(nwb_copy_file_name, mode='w', manager=nwb_io.manager) as io:
                io.write(nwb_file)
                io.close()
                nwb_io.close()
                return nwb_copy_file_name
        else:
            if verbose: print("Writing NWB file and closing.")
            nwb_io.write(nwb_file)
            nwb_io.close()
            return nwb_file_name
    else:
        if verbose: print("Returning NWB file variable")
        return nwb_file



def NWBWaveletTransform(nwb_file, signal_name=None, freq=None, width=None, elec_ids=None, wavelet_type='morlet', downsample=None, copy_file=False, verbose=True):
    """
    Does a wavelet transform on electrophysiology data from NWB file. Saves in a new processing module
    in the same NWB file.

    :param nwb_file: {NWB file or str} Either a handle to the open NWB file or a string with the file name.
    :param signal_name: {str} Name of the signal to process. This should be the name of the ElectricalSeries.  Will search for this in raw or LFP module.
    :param freq: [optional] {list of floats} List of center frequencies for the wavelets.  Default is 1 to 100 Hz, every 1 Hz.
    :param width: [optional] {list of floats} Width of the wavelet to be used (measured in # of cycles in time). Default is 2.5.
    :param elec_ids: [optional] {List of integers} List of electrodes to decompose. Default is all electrodes.
    :param wavelet_type: [optional] {str} Type of wavelet to use. Currently only 'morlet' is valid.
    :param downsample: [optional] {int} If specified, downsample signal by this ratio. Default is None (keep original)
    :param verbose: [optional] {bool} Whether to print updates of progress. Default is False.
    :return: {NWB file or str} Depending on type of nwb_file input passed, either returns an open NWB file handle or the string of the file name.
    """

    # Check to see if user passed a filename or file
    if not nwb_file:
        raise ValueError("Must pass a NWB file.")
    if type(nwb_file) == file.NWBFile:
        if verbose: print("Using passed NWB file.")
    elif type(nwb_file) == str:
        nwb_file_name = nwb_file
        [nwb_file, nwb_io] = OpenNWBFile(nwb_file_name, verbose=verbose)
    else:
        raise ValueError("NWB file was not valid.")

    # Check to see if we should do a copy of the file first
    if copy_file:
        nwb_file = nwb_file.copy()

    # Parse name of signal we want to transform, grab electrical series
    if signal_name in nwb_file.acquisition.keys():
        # Then use raw data
        cur_es = nwb_file.acquisition[signal_name]
        cur_es_data = nwb_file.acquisition[signal_name].data
        es_elec_ids = nwb_file.acquisition[signal_name].electrodes.table[:,0]
        cur_es_rate = nwb_file.acquisition[signal_name].rate
        cur_es_starting_time = nwb_file.acquisition[signal_name].starting_time
    elif ('ecephys' in nwb_file.processing.keys()) and (signal_name in nwb_file.processing['ecephys'].data_interfaces.keys()):
        # Then use LFP data from this interface -- concatenate them all
        cur_es = nwb_file.processing['ecephys'][signal_name].electrical_series
        cur_es_data = []
        for cur_key in cur_es.keys():
            if not cur_es_data:
                cur_es_data = cur_es[cur_key].data
                es_elec_ids = cur_es[cur_key].electrodes.table[:,0]
                cur_es_rate = cur_es[cur_key].rate
                cur_es_starting_time = cur_es[cur_key].starting_time
            else:
                cur_es_data = np.concatenate((cur_es_data, cur_es[cur_key].data), axis=1)
                es_elec_ids = np.concatenate((es_elec_ids, cur_es[cur_key].electrodes.table[:,0]))

    # Parse wavelet type
    if not wavelet_type:
        wavelet_type = 'morlet'
    wavelet_type = wavelet_type.lower()

    # Parse inputs on wavelets
    if not freq:
        freq = np.arange(1, 101, 1)
    if not width:
        width = 6

    # Convert width (# of cycles) to standard deviation
    if wavelet_type == 'morlet':
        wavelet_widths = width * cur_es_rate / (2 * freq * np.pi)

    # Plot all of the wavelets to a figure
    if verbose:
        fig = plt.figure()
        max_width = np.round(np.max(wavelet_widths)*width*1.5/2)
        time = np.arange(-max_width, max_width, 1)/cur_es_rate
        for cur_f_ind, cur_f in enumerate(freq):
            plt.plot(time, signal.morlet2(2*max_width, wavelet_widths[cur_f_ind], w=width), label="%d Hz" % (cur_f))
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude (arb)')
        plt.title('Time series of wavelets')
        fig.axes[0].legend()
        plt.show()

    # Create new processing module for transformed data
    # Make sure extracellular ephys module exists in the file
    if 'ecephys' not in nwb_file.processing.keys():
        if verbose: print("Electrophysiology processing module does not exist. Creating.")
        signal_module = ProcessingModule(name='ecephys',
                                         description="Processing module for electrophysiology signals.")
        nwb_file.add_processing_module(signal_module)

    # Check our list of electrode ids
    if type(elec_ids) is type(None):
        elec_ids = es_elec_ids
    if type(elec_ids) is not np.ndarray:
        elec_ids = np.array(elec_ids)
    elec_ids = elec_ids[np.in1d(elec_ids, es_elec_ids)]

    # Initialize the CWT data matrix
    num_samples = cur_es_data.shape[0]
    cur_rate = cur_es_rate
    if np.isreal(downsample):
        downsample = round(float(downsample))
        if downsample > 1:
            num_samples = int(round(cur_es_data.shape[0]/downsample))
            cur_rate = cur_rate/downsample
    cwt_data = np.NaN * np.ones((num_samples, len(elec_ids), len(freq)), dtype=np.complex128)

    # Apply wavelets to each electrode in turn
    if verbose: print("Wavelet transforming signal (wavelet=%s):" % (wavelet_type))
    for cur_elec in elec_ids:
        cur_sig_ind = np.where(es_elec_ids==cur_elec)[0][0]
        if verbose: print("\tTransforming signal from electrode %d (signal %d)..." % (cur_elec, cur_sig_ind), end='')
        # Apply wavelet to ephys signal
        if wavelet_type == 'morlet':
            cur_cwt = signal.cwt(cur_es_data[:, cur_sig_ind], signal.morlet2, wavelet_widths, w=width)
        else:
            raise ValueError("Requested wavelet type is not valid.")
        # Downsample to new frequency
        if (downsample > 1):
            cur_cwt = cur_cwt[:,0::downsample]
        # Add to overall data
        cwt_data[:, cur_sig_ind, :] = cur_cwt.transpose()
        if verbose: print("done.")

    # Create spectral decomposition data interface
    cwt_decomp = DecompositionSeries('CWT of %s' % (signal_name),
                                     cwt_data,
                                     metric='complex',
                                     description="Wavelet decomposition (%s, width=%d) of ephys signal from %s. Electrodes: [%s]" %
                                                 (wavelet_type, width, signal_name, elec_ids),
                                     unit='Arbitrary (complex)',
                                     source_timeseries=cur_es,
                                     starting_time=cur_es.starting_time,
                                     rate=cur_rate,
                                     comments='This is the complex transform. Further analysis needed for getting magnitude of phase')
    for cur_f in freq:
        cwt_decomp.add_band(band_name="%d Hz",
                            band_limits=[float(cur_f) - 3*float(cur_f)/float(width), float(cur_f) + 3*float(cur_f)/float(width)],
                            band_mean=float(cur_f),
                            band_stdev=float(cur_f)/float(width))

    # Add data interface to behavior module in NWB file
    if verbose: print("Adding wavelet-transformed data interface to electrophysiology module.")
    nwb_file.processing['ecephys'].add(cwt_decomp)

    # Write the file
    if nwb_io is not None:
        if copy_file:
            nwb_copy_file_name = path.splitext(nwb_file_name)[0] + '_CWTProcessed' + path.splitext(nwb_file_name)[1]
            if verbose: print("Did a shallow copy of NWB file, so saving to new file: %s" % (nwb_copy_file_name))
            with NWBHDF5IO(nwb_copy_file_name, mode='w', manager=nwb_io.manager) as io:
                io.write(nwb_file)
                io.close()
                return nwb_copy_file_name
        else:
            if verbose: print("Writing NWB file and closing.")
            nwb_io.write(nwb_file)
            nwb_io.close()
            return nwb_file_name
    else:
        if verbose: print("Returning NWB file variable")
        return nwb_file


