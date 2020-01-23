# Functions for importing blackrock files into NWB format
from qtpy.QtWidgets import QFileDialog, QApplication
from os import getcwd, path, listdir
import sys
from datetime import datetime
from dateutil.tz import tzlocal
import numpy as np
import re

from pynwb import NWBFile, NWBHDF5IO, ProcessingModule, file
from pynwb.ecephys import ElectricalSeries, FilteredEphys, SpikeEventSeries, EventWaveform
from pynwb.behavior import BehavioralTimeSeries, EyeTracking, PupilTracking
from pynwb.misc import AbstractFeatureSeries

import tdt


def InitializeNWBFromTDT(tdt_block_dir, nwb_file_name='', experiment=None, subject=None,
                               electrode_info=None, electrode_group=None, notes=None, overwrite=False, verbose=False):
    """
        Initializes a NWB file to copy over data from Blackrock file.  User specified information is
        added to NWB file as metadata about data.

        :param tdt_block_dir: {str} full path of TDT files to convert to NWB. If empty, will open dialog.
        :param nwb_file_name: [optional] {str} full path of NWB file to export to. Default is to change TDT extension to 'nwb'
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :param experiment: [optional] {dict} Dictionary that contains information about experimenter, name, animal, lab, and institution
        :param subject: [optional] {dict} Dictionary that contains information about subject
        :param electrode_info: [optional] {list} List of dictionaries carrying information about electrodes
        :param electrode_group: [optional] {list} List of dictionaries carrying information about electrode groupings
        :param notes: [optional] {str} Notes relevant for the dataset
        :param overwrite: [optional] {bool} Whether to overwrite the NWB file if it already exists
        :return: {str} filename of NWB file (empty if error)
        """
    # Check to see if user specified a TDT filename
    if not tdt_block_dir:  # no file name passed
        # Ask user to specify a file
        if 'app' not in locals():
            app = QApplication([])
        tdt_block_dir = QFileDialog.getExistingDirectory(QFileDialog(), 'Select Directory', getcwd())

    # Check to see if valid nwb_file_name is passed
    tdt_tsq_files = [f for f in listdir(tdt_block_dir) if path.splitext(f)[1] == '.tsq']
    if not nwb_file_name:
        nwb_file_name = tdt_block_dir + path.splitext(tdt_tsq_files[0])[0] + '.nwb'
    if verbose: print("Writing to NWB data file %s" % (nwb_file_name))

    # Initialize the TDT file, get header and basic information
    try:
        # Requires the raw data to be imported
        tdt_header = tdt.read_block(tdt_block_dir, headers=1)
        tdt_data = tdt.read_block(tdt_block_dir, evtype=['epocs'])
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        raise FileNotFoundError("Couldn't open TDT file. Error: %s" % e)

    # Process experiment information from inputs
    if (not experiment) or not (type(experiment) is dict):
        if verbose: print("Using default experimental details.")
        nwb_session_description = 'Experimental details not provided'
        experimenter = ''
        lab = ''
        institution = ''
    else:
        if 'name' not in experiment.keys():
            experiment['name'] = ''
        if 'animal' not in experiment.keys():
            experiment['animal'] = ''
        if 'experimenter' not in experiment.keys():
            experiment['experimenter'] = ''
        if 'lab' not in experiment.keys():
            experiment['lab'] = ''
        if 'institution' not in experiment.keys():
            experiment['institution'] = ''
        experimenter = experiment['experimenter']
        lab = experiment['lab']
        institution = experiment['institution']
        nwb_session_description = 'Experiment: ' + experiment['name'] + \
                                  ', Animal: ' + experiment['animal'] + \
                                  ', Date: ' + tdt_data.info.start_date.strftime('%m/%d/%Y')  # description of the recording session
        if verbose: print("Experiment description: %s" % (nwb_session_description))

    # Process subject information
    if (not subject) or not (type(subject) is dict):
        cur_subject = None
    else:
        try:
            # Create a subject object from the passed dictionary
            cur_subject = file.Subject(**subject)
        except:  # catch *all* exceptions
            e = sys.exc_info()[0]
            raise ValueError("Couldn't create subject object. Error: %s" % (e))

    # Define the NWB file identifier as the filename, as this should be unique
    nwb_identifier = path.split(tdt_tsq_files[0])[1]  # unique file identifier

    # Initialize the NWB file
    nwb_io = []
    try:
        if not overwrite and path.isfile(nwb_file_name):
            raise FileExistsError(
                "Can't update basic attributes of NWB file. Create new NWB file if you need to change attributes.")
        else:
            if verbose: print("Creating new NWB file %s..." % (nwb_file_name), end='')
            nwb_file_append = False
            nwb_io = NWBHDF5IO(nwb_file_name, mode='w')
            tdt_notes = ''
            if 'Note' in tdt_data.epocs.keys():
                for cur_note in tdt_data.epocs.Note:
                    tdt_notes = tdt_notes + cur_note.notes + '; '
            nwb_file = NWBFile(nwb_session_description,
                               nwb_identifier,
                               tdt_data.info.start_date,
                               file_create_date=datetime.now(tzlocal()),
                               experimenter=experimenter,
                               subject=cur_subject,
                               lab=lab,
                               institution=institution,
                               source_script='ConvertTDTToNWB.py',
                               source_script_file_name='ConvertTDTToNWB.py',
                               notes=notes + tdt_notes)
            if verbose: print("done.")
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        if nwb_io: nwb_io.close()
        tdt_data.clear()
        raise FileExistsError("Couldn't open NWB file. Error: %s" % e)
    if verbose: print("Created NWB file.")

    # Create device in the NWB file
    device = nwb_file.create_device(name='TDT')
    if verbose: print("Created device.")

    # Make sure electrode input is provided and valid
    if (not electrode_info) or not (type(electrode_info) is list):
        if verbose: print("Creating electrode information from directory.")

        # Grab information about SEV files
        tdt_sev_files = [f for f in listdir(tdt_block_dir) if path.splitext(f)[1] == '.sev']
        electrode_list = []
        for sev_ind, sev_file in enumerate(tdt_sev_files):
            sev_match = re.search(".+_RAW(\\d+)_Ch(\\d+).sev", sev_file)
            electrode_list.append((int(sev_match.groups()[0])-1)*16 + int(sev_match.groups()[1]))
        electrode_list = np.sort(electrode_list)

        # Create electrode information for each electrode
        electrode_info = []
        for cur_elec in electrode_list:
            electrode_info.append({'x': 0.0, 'y': 0.0, 'z': 0.0, 'impedance': -1.0,
                                   'location': 'unknown', 'group': 'electrodes',
                                   'id': cur_elec,
                                   'filtering': 'Unknown'})
        if verbose: print("\tCreated %d electrodes." % (len(electrode_info)))

    # Make sure electrode group input is provided and valid
    if verbose: print("Processing electrode groups.")
    default_electrode_group = {'name': 'default electrode group',
                               'description': 'Generic electrode group for ungrouped electrodes.',
                               'location': 'Unknown'}
    if (not electrode_group) or not (type(electrode_group) is list):
        electrode_group = [default_electrode_group]
    else:
        electrode_group.insert(0, default_electrode_group)

    # Fill in any missing information about electrodes
    ovr_elec_group_list = []  # will be used to keep track of electrode group list
    if verbose: print("Processing electrode information.")
    for cur_elec_info_ind, cur_elec_info in enumerate(electrode_info):
        # If id is not provided, then throw an error
        if ('id' not in cur_elec_info.keys()) or (not cur_elec_info['id']):
            tdt_data.clear()
            nwb_io.close()
            raise ValueError("Couldn't process electrode inputs. Must provide IDs for all electrodes or pass empty array.")
        # If filtering is not provided, add it from the file (can't with TDT)
        if ('filtering' not in cur_elec_info.keys()) or (not cur_elec_info['filtering']):
            electrode_info[cur_elec_info_ind]['filtering'] = 'Unknown'
        # If other variables are not provided, just initialize them to defaults
        if ('x' not in cur_elec_info.keys()) or (not cur_elec_info['x']):
            electrode_info[cur_elec_info_ind]['x'] = 0.0
        if ('y' not in cur_elec_info.keys()) or (not cur_elec_info['y']):
            electrode_info[cur_elec_info_ind]['y'] = 0.0
        if ('z' not in cur_elec_info.keys()) or (not cur_elec_info['z']):
            electrode_info[cur_elec_info_ind]['z'] = 0.0
        if ('location' not in cur_elec_info.keys()) or (not cur_elec_info['location']):
            electrode_info[cur_elec_info_ind]['location'] = 'unknown'
        if ('impedance' not in cur_elec_info.keys()) or (not cur_elec_info['impedance']):
            electrode_info[cur_elec_info_ind]['impedance'] = float(-1.0)
        if ('group' not in cur_elec_info.keys()) or (not cur_elec_info['group']):
            electrode_info[cur_elec_info_ind]['group'] = 'default electrode group'
        # Find the index to the electrode group
        grp_elec_ind = 0
        while (grp_elec_ind < len(electrode_group)) and (
                cur_elec_info['group'] != electrode_group[grp_elec_ind]['name']):
            grp_elec_ind = grp_elec_ind + 1
        # If we made it past the end, this electrode group doesn't exist, put it in the default
        if grp_elec_ind >= len(electrode_group):
            grp_elec_ind = -1
        # Save the index into group dictionary for easy reference
        electrode_info[cur_elec_info_ind]['electrode_group_index'] = grp_elec_ind
        ovr_elec_group_list.append(grp_elec_ind)

    # Create electrode groups
    nwb_electrode_group = []
    if verbose: print("Creating electrode groups.")
    for cur_elec_group_ind, cur_elec_group in enumerate(electrode_group):
        # Create and add to our list
        nwb_electrode_group.append(nwb_file.create_electrode_group(cur_elec_group['name'],
                                                                   description=cur_elec_group['description'],
                                                                   location=cur_elec_group['location'], device=device))
    if verbose: print("\tCreated %d electrode groups." % (len(electrode_group)))

    # Create electrodes in NWB file
    if verbose: print("Adding electrodes to NWB.")
    nwb_file_elec_list = []
    for cur_elec_ind, cur_elec in enumerate(electrode_info):
        # Add electrode to NWB file
        nwb_file.add_electrode(id=cur_elec['id'], x=cur_elec['x'], y=cur_elec['y'], z=cur_elec['z'],
                               imp=cur_elec['impedance'], location=cur_elec['location'],
                               filtering=cur_elec['filtering'],
                               group=nwb_electrode_group[cur_elec['electrode_group_index']])
        # Keep track of electrodes entered so we can index for the electrode table
        nwb_file_elec_list.append(cur_elec['id'])
    if verbose: print("\tAdded %d electrodes." % (len(electrode_info)))

    # Close NSx file
    tdt_data.clear()

    # Write to file and close
    if verbose: print("\tWriting NWB file and closing.")
    nwb_io.write(nwb_file)
    nwb_io.close()

    return nwb_file_name


def AddTDTRawDataToNWB(tdt_block_dir, nwb_file_name='', verbose=False, elec_ids=None):
    """
        Copies raw electrophysiological data from TDT format to Neurodata Without Borders (NWB) format.
        Typically, the NWB file will already be initialized by calling InitializeNWBFromTDT. Multiple
        electrodes can be added at once or with separate calls to this function.

        :param tdt_block_dir: {str} full path of TDT files to convert to NWB. If empty, will open dialog.
        :param nwb_file_name: [optional] {str} full path of NWB file to export to. Default is to change blackrock extension to 'nwb'
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :param elec_ids: [optional] {list} List of electrode IDs to copy over. If empty, all are copied
        :return: {str} filename of NWB file (empty if error)
        """
    # Check to see if user specified a TDT filename
    if not tdt_block_dir:  # no file name passed
        # Ask user to specify a file
        if 'app' not in locals():
            app = QApplication([])
        tdt_block_dir = QFileDialog.getExistingDirectory(QFileDialog(), 'Select Directory', getcwd())

    # Check to see if valid nwb_file_name is passed
    tdt_tsq_files = [f for f in listdir(tdt_block_dir) if path.splitext(f)[1] == '.tsq']
    if not nwb_file_name:
        nwb_file_name = path.splitext(tdt_tsq_files[0])[0] + '.nwb'
    if verbose: print("Writing to NWB data file %s" % (nwb_file_name))

    # Initialize the TDT file
    try:
        # Requires the raw data to be imported
        tdt_header = tdt.read_block(tdt_block_dir, headers=1)
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        raise FileNotFoundError("Couldn't open TDT file. Error: {:s}".format(e))

    # Initialize the NWB file
    nwb_io = []
    try:
        if not path.isfile(nwb_file_name):
            # Initialize NWB file
            if verbose: print("NWB file doesn't exist. Creating new one: %s..." % (nwb_file_name))
            InitializeNWBFromTDT(tdt_block_dir, nwb_file_name, verbose=verbose)

        # Append to existing file
        if verbose: print("Opening NWB file %s..." % (nwb_file_name), end='')
        nwb_file_append = True
        nwb_io = NWBHDF5IO(nwb_file_name, mode='a')
        nwb_file = nwb_io.read()
        if verbose: print("done.")
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        if nwb_io: nwb_io.close()
        raise FileExistsError("Couldn't open NWB file. Error: %s" % e)

    # Validate the elec_ids list
    if not elec_ids:
        # Grab electroe list from NWB file
        elec_ids = nwb_file.electrodes.id[:]

    # Load data for all of the electrodes
    if verbose: print("Loading raw data from TDT file...")
    chan_data = []
    starting_time = 0
    fs = 24414.0625
    for cur_elec_id in elec_ids:
        # Convert electrode id into RAW and channel
        cur_raw = np.floor((cur_elec_id-1)/16) + 1
        cur_raw = 'RAW{:1.0f}'.format(cur_raw)
        cur_chan = ((cur_elec_id-1) % 16) + 1
        print('{:s}_{:1.0f}'.format(cur_raw, cur_chan))

        #Read this channel data in from file
        cur_chan_data = tdt.read_sev(tdt_block_dir, channel=cur_chan, event_name=cur_raw)
        fs = cur_chan_data[cur_raw].fs
        starting_time = cur_chan_data.time_ranges[0][0]
        cur_chan_data = cur_chan_data[cur_raw].data
        cur_chan_data = np.reshape(cur_chan_data, [len(cur_chan_data), 1])
        if len(chan_data) == 0:
            chan_data = cur_chan_data
        else:
            chan_data = np.concatenate([chan_data, cur_chan_data], axis=1)

        if verbose: print("\tLoaded electrode %d..." % (cur_elec_id))

    # Get resolution of signal in V
    if verbose: print("Processing electrode information.")
    V_resolution = 0.0 #TODO: See if we can get resolution information from TDT?

    # Create electrode table for new electrodes
    elec_table_ind = np.ones((len(elec_ids), 1)) * np.NaN
    for cur_elec_ind, cur_elec in enumerate(elec_ids):
        # In order to create electrode table, we have to have indexes for each channel into electrode list in NWB file
        cur_elec_table_ind = 0
        while (cur_elec_table_ind < len(nwb_file.electrodes)) and (
                nwb_file.electrodes[cur_elec_table_ind, 0] != cur_elec):
            cur_elec_table_ind = cur_elec_table_ind + 1
        if cur_elec_table_ind >= len(nwb_file.electrodes):
            raise ValueError("Couldn't find electrode %d in NWB file list." % (cur_elec))
        elec_table_ind[cur_elec_ind] = cur_elec_table_ind
    elec_table_ind = elec_table_ind.transpose().tolist()[0]  # Convert to list for create_electrode_table_region
    electrode_table_region = nwb_file.create_electrode_table_region(elec_table_ind, "Electrodes %s" % (elec_ids))

    # Write raw data for the electrode
    if verbose: print("\tAdding to NWB...", end='')
    ephys_ts = ElectricalSeries('Raw Data for Channels %s' % (elec_ids),
                                chan_data,
                                electrode_table_region,
                                starting_time=starting_time,
                                rate=fs,
                                resolution=V_resolution,
                                conversion=float(1.0), #TODO: Check that TDT does read in V, as expected by NWB
                                comments="",
                                description="Raw data from TDT file.")
    nwb_file.add_acquisition(ephys_ts)
    if verbose: print("done")

    # Write the file
    if verbose: print("\tWriting NWB file and closing.")
    nwb_io.write(nwb_file)
    nwb_io.close()

    return nwb_file_name

def AddTDTAnalogDataToNWB(tdt_block_dir, nwb_file_name='', signal_info=None, module_name='behavior', verbose=False):
    """
        Copies analog (continuous) data from the specified Blackrock file to Neurodata Without Borders (NWB) file.
        This is usually continuous signals about behavior (joystick position, screen refreshes, etc).
        User should provide information about the signals to help end users understand what is in each signal.
        Multiple calls can be used to load multiple data signals from multiple files.

        Typically, the NWB file will already be initialized by calling InitializeNWBFromBlackrock.

        :param tdt_block_dir: {str} full path of TDT files to convert to NWB. If empty, will open dialog.
        :param nwb_file_name: [optional] {str} full path of NWB file to export to. Default is to change blackrock extension to 'nwb'
        :param signal_info: [optional] {list} List of dictionaries with information about the signals to save.
        :param module_name: [optional] {str} Name of module to store data.  Usually 'behavior' but could also be 'ecephys' or 'misc'
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :return: {str} filename of NWB file (empty if error)
        """
    # Check to see if user specified a TDT filename
    if not tdt_block_dir:  # no file name passed
        # Ask user to specify a file
        if 'app' not in locals():
            app = QApplication([])
        tdt_block_dir = QFileDialog.getExistingDirectory(QFileDialog(), 'Select Directory', getcwd())

    # Check to see if valid nwb_file_name is passed
    tdt_tsq_files = [f for f in listdir(tdt_block_dir) if path.splitext(f)[1] == '.tsq']
    if not nwb_file_name:
        nwb_file_name = path.splitext(tdt_tsq_files[0])[0] + '.nwb'
    if verbose: print("Writing to NWB data file %s" % (nwb_file_name))

    # Initialize the TDT file
    try:
        # Requires the raw data to be imported
        tdt_header = tdt.read_block(tdt_block_dir, headers=1)
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        raise FileNotFoundError("Couldn't open TDT file. Error: {:s}".format(e))

    # Initialize the NWB file
    nwb_io = []
    try:
        if not path.isfile(nwb_file_name):
            # Initialize NWB file
            if verbose: print("NWB file doesn't exist. Creating new one: %s..." % (nwb_file_name))
            InitializeNWBFromTDT(tdt_block_dir, nwb_file_name, verbose=verbose)

        # Append to existing file
        if verbose: print("Opening NWB file %s..." % (nwb_file_name), end='')
        nwb_file_append = True
        nwb_io = NWBHDF5IO(nwb_file_name, mode='a')
        nwb_file = nwb_io.read()
        if verbose: print("done.")
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        if nwb_io: nwb_io.close()
        raise FileExistsError("Couldn't open NWB file. Error: %s" % e)

    # Make sure module name is either behavior or misc
    module_name = module_name.lower()
    if (module_name != 'behavior') and (module_name != 'misc'):
        raise ValueError("Module type must either be 'behavior' or 'misc'.")

    # Parse the signal_info list
    if not signal_info:
        raise ValueError("Must specify signals to load.")
    elec_ids = []
    for cur_signal_ind, cur_signal_info in enumerate(signal_info):
        if 'label' not in cur_signal_info.keys():
            raise ValueError("Signal information must have a label for each signal.")
        if 'name' not in cur_signal_info.keys():
            raise ValueError("Signal information must have a name for each signal. (Should be user-understandable)")
        if 'comments' not in cur_signal_info.keys():
            signal_info[cur_signal_ind]['comments'] = ''
        # Find electrode IDs for this signal
        if ('elec_id' not in cur_signal_info.keys()) or (not cur_signal_info['elec_id']):
            # Loop through and grab all signals of type 'streams' that aren't RAW data
            signal_info[cur_signal_ind]['elec_id'] = []
            for cur_store in tdt_header.stores.keys():
                # Grab all 'streams' but ignore RAWs
                if (tdt_header.stores[cur_store].type_str == 'streams') and (cur_store[0:3] != 'RAW'):
                    signal_info[cur_signal_ind]['elec_id'].append(cur_store)

    # Create processing module for saving data
    if module_name not in nwb_file.processing.keys():
        if verbose: print("Specified processing module (%s) does not exist. Creating." % (module_name))
        signal_module = ProcessingModule(name=module_name,
                                         description="Processing module for continuous signal data from %s." % (
                                             path.split(tdt_tsq_files[0])[1]))
        nwb_file.add_processing_module(signal_module)

    # Create data interface for the analog signals
    signal_info_str = signal_info[0]['name']
    for i in range(1, len(signal_info)):
        signal_info_str = signal_info_str + ", " + signal_info[i]['name']
    if verbose: print("Creating %s data interface for signals %s." % (module_name, signal_info_str))
    if module_name == 'behavior':
        cur_data_interface = BehavioralTimeSeries(name="Analog signals (" + signal_info_str + ")")
    elif module_name == 'misc':
        cur_data_interface = AbstractFeatureSeries(name="Analog signals (" + signal_info_str + ")")
    else:
        raise ValueError("Module type must either be 'behavior' or 'misc'.")

    if verbose: print("Adding signals...")
    for cur_signal_ind, cur_signal_info in enumerate(signal_info):
        # Get data from file
        analog_data = []
        analog_fs = []
        analog_start_time = []
        for cur_elec_id in cur_signal_info['elec_id']:
            cur_data = tdt.read_block(tdt_block_dir, store=cur_elec_id)
            cur_analog_data = cur_data.streams[cur_elec_id].data
            cur_analog_data = np.reshape(cur_analog_data, [len(cur_analog_data), 1])
            if len(analog_data) == 0:
                analog_fs = cur_data.streams[cur_elec_id].fs
                analog_start_time = cur_data.streams[cur_elec_id].start_time
                analog_data = cur_analog_data
            else:
                analog_fs.append(cur_data.streams[cur_elec_id].fs)
                analog_start_time.append(cur_data.streams[cur_elec_id].start_time)
                analog_data = np.concatenate([analog_data, cur_analog_data], axis=1)

        # Make sure all of the fs and start_times are the same
        analog_start_time = np.unique(analog_start_time)
        analog_fs = np.unique(analog_fs)
        if len(analog_start_time) != 1 or len(analog_fs) != 1:
            raise ValueError('Start time and sampling frequency need to be the same for signals to be combined.')

        # Create time series
        cur_data_interface.create_timeseries(name=cur_signal_info['name'],
                                             data=analog_data,
                                             comments=cur_signal_info['comments'],
                                             unit="V", #TODO: Check that this is correct for TDT
                                             resolution=1.0, #TODO: Can we get this from TDT?
                                             conversion=0.001, #TODO: Check what the correct conversion is for TDT
                                             starting_time=analog_start_time[0],
                                             rate=analog_fs[0],
                                             description="Signal %s from %s." % (
                                                 cur_signal_info['label'],
                                                 path.split(tdt_tsq_files[0])[1]))
        if verbose: print("\tAdded %s." % (cur_signal_info['label']))

    # Add data interface to module in NWB file
    if verbose: print("Adding data interface to module.")
    nwb_file.processing[module_name].add(cur_data_interface)

    # Write the file
    if verbose: print("Writing NWB file and closing.")
    nwb_io.write(nwb_file)
    nwb_io.close()

    return nwb_file_name