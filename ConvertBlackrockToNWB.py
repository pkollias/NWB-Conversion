# Functions for importing blackrock files into NWB format
from qtpy.QtWidgets import QFileDialog, QApplication
from os import getcwd, path
import warnings
import sys
from datetime import datetime
from dateutil.tz import tzlocal
import numpy as np

from pynwb import NWBFile, NWBHDF5IO, ProcessingModule, file
from pynwb.ecephys import ElectricalSeries, FilteredEphys, SpikeEventSeries, EventWaveform
from pynwb.behavior import BehavioralTimeSeries, EyeTracking, PupilTracking
from pynwb.misc import AbstractFeatureSeries

import brpylib


def InitializeNWBFromBlackrock(blk_file_name, nwb_file_name='', experiment=None, subject=None,
                               electrode_info=None, electrode_group=None, notes=None, overwrite=False, verbose=False):
    """
        Initializes a NWB file to copy over data from Blackrock file.  User specified information is
        added to NWB file as metadata about data.

        :param blk_file_name: {str} full path of Blackrock file to convert to NWB. If empty, will open dialog.
        :param nwb_file_name: [optional] {str} full path of NWB file to export to. Default is to change blackrock extension to 'nwb'
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :param experiment: [optional] {dict} Dictionary that contains information about experimenter, name, animal, lab, and institution
        :param subject: [optional] {dict} Dictionary that contains information about subject
        :param electrode_info: [optional] {list} List of dictionaries carrying information about electrodes
        :param electrode_group: [optional] {list} List of dictionaries carrying information about electrode groupings
        :param notes: [optional] {str} Notes relevant for the dataset
        :param overwrite: [optional] {bool} Whether to overwrite the NWB file if it already exists
        :return: {str} filename of NWB file (empty if error)
        """
    # Check to see if user specified a Blackrock filename
    if not blk_file_name:  # no file name passed
        # Ask user to specify a file
        if 'app' not in locals():
            app = QApplication([])
        blk_file_name = QFileDialog.getOpenFileName(QFileDialog(), "Select File", getcwd(),
                                                    "Blackrock Data File (*.ns6)")
        blk_file_name = blk_file_name[0]

    # Check to see if valid nwb_file_name is passed
    if not nwb_file_name:
        nwb_file_name = path.splitext(blk_file_name)[0] + '.nwb'
    if verbose: print("Writing to NWB data file: %s" % (nwb_file_name))

    # Initialize the Blackrock file
    try:
        # Requires the raw data to be imported
        if verbose: print("Initializing read from Blackrock file: %s" % (blk_file_name))
        nsx_file = brpylib.NsxFile(blk_file_name)
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        raise FileNotFoundError("Couldn't open Blackrock file. Error: %s" % e)

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
                                  ', Date: ' + nsx_file.basic_header['TimeOrigin'].strftime(
            '%m/%d/%Y')  # description of the recording session
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
    nwb_identifier = path.split(blk_file_name)[1]  # unique file identifier

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
            nwb_file = NWBFile(nwb_session_description, nwb_identifier, nsx_file.basic_header['TimeOrigin'],
                               file_create_date=datetime.now(tzlocal()),
                               experimenter=experimenter,
                               subject=cur_subject,
                               lab=lab,
                               institution=institution,
                               source_script='ConvertBlackrockToNWB.py',
                               source_script_file_name='ConvertBlackrockToNWB.py',
                               notes=notes + "Comments from Blackrock: " + nsx_file.basic_header['Comment'])
            if verbose: print("done.")
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        if nwb_io: nwb_io.close()
        nsx_file.datafile.close()
        raise FileExistsError("Couldn't open NWB file. Error: %s" % e)
    if verbose: print("Created NWB file.")

    # Create device in the NWB file
    device = nwb_file.create_device(
        name='Blackrock (' + nsx_file.basic_header['FileTypeID'] + ', ' + nsx_file.basic_header['FileSpec'] + ')')
    if verbose: print("Created device.")

    # Make sure electrode input is provided and valid
    if (not electrode_info) or not (type(electrode_info) is list):
        if verbose: print("Creating electrode information from file.")
        # Create electrode information for each electrode
        electrode_info = []
        for cur_elec in range(nsx_file.basic_header['ChannelCount']):
            electrode_info.append({'x': 0.0, 'y': 0.0, 'z': 0.0, 'impedance': -1.0,
                                   'location': 'unknown', 'group': 'electrodes',
                                   'id': nsx_file.extended_headers[cur_elec]['ElectrodeID'],
                                   'filtering': 'High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)' %
                                                (nsx_file.extended_headers[0]['HighFreqCorner'],
                                                 nsx_file.extended_headers[0]['HighFreqOrder'],
                                                 nsx_file.extended_headers[0]['HighFreqType'],
                                                 nsx_file.extended_headers[0]['LowFreqCorner'],
                                                 nsx_file.extended_headers[0]['LowFreqOrder'],
                                                 nsx_file.extended_headers[0]['LowFreqType'])})
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
            nsx_file.datafile.close()
            nwb_io.close()
            raise ValueError(
                "Couldn't process electrode inputs. Must provide IDs for all electrodes or pass empty array.")
        # If filtering is not provided, add it from NSX file
        if ('filtering' not in cur_elec_info.keys()) or (not cur_elec_info['filtering']):
            electrode_info[cur_elec_info_ind][
                'filtering'] = 'High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)' % (
                nsx_file.extended_headers[0]['HighFreqCorner'], nsx_file.extended_headers[0]['HighFreqOrder'],
                nsx_file.extended_headers[0]['HighFreqType'], nsx_file.extended_headers[0]['LowFreqCorner'],
                nsx_file.extended_headers[0]['LowFreqOrder'], nsx_file.extended_headers[0]['LowFreqType'])
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
        if grp_elec_ind >= len(electrode_group): grp_elec_ind = -1
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
    nsx_file.datafile.close()

    # Write to file and close
    if verbose: print("\tWriting NWB file and closing.")
    nwb_io.write(nwb_file)
    nwb_io.close()

    return nwb_file_name


def AddBlackrockRawDataToNWB(blk_file_name, nwb_file_name='', verbose=False, elec_ids=None):
    """
        Copies raw electrophysiological data from Blackrock format to Neurodata Without Borders (NWB) format.
        Typically, the NWB file will already be initialized by calling InitializeNWBFromBlackrock. Multiple
        electrodes can be added at once or with separate calls to this function.

        :param blk_file_name: {str} full path of Blackrock file to convert to NWB. If empty, will open dialog.
        :param nwb_file_name: [optional] {str} full path of NWB file to export to. Default is to change blackrock extension to 'nwb'
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :param elec_ids: [optional] {list} List of electrode IDs to copy over. If empty, all are copied
        :return: {str} filename of NWB file (empty if error)
        """
    # Check to see if user specified a Blackrock filename
    if not blk_file_name:  # no file name passed
        # Ask user to specify a file
        if 'app' not in locals():
            app = QApplication([])
        blk_file_name = QFileDialog.getOpenFileName(QFileDialog(), "Select File", getcwd(),
                                                    "Blackrock Data File (*.ns6)")
        blk_file_name = blk_file_name[0]

    # Check to see if valid nwb_file_name is passed
    if not nwb_file_name:
        nwb_file_name = path.splitext(blk_file_name)[0] + '.nwb'
    if verbose: print("Writing to NWB data file %s" % (nwb_file_name))

    # Initialize the Blackrock file
    try:
        # Requires the raw data to be imported
        nsx_file = brpylib.NsxFile(blk_file_name)
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        raise FileNotFoundError("Couldn't open Blackrock file. Error: %s" % e)

    # Initialize the NWB file
    nwb_io = []
    try:
        if not path.isfile(nwb_file_name):
            # Initialize NWB file
            if verbose: print("NWB file doesn't exist. Creating new one: %s..." % (nwb_file_name))
            InitializeNWBFromBlackrock(blk_file_name, nwb_file_name, verbose=True)

        # Append to existing file
        if verbose: print("Opening NWB file %s..." % (nwb_file_name), end='')
        nwb_file_append = True
        nwb_io = NWBHDF5IO(nwb_file_name, mode='a')
        nwb_file = nwb_io.read()
        if verbose: print("done.")
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        if nwb_io: nwb_io.close()
        nsx_file.datafile.close()
        raise FileExistsError("Couldn't open NWB file. Error: %s" % e)

    # Validate the elec_ids list
    if not elec_ids:
        # Grab all of the data
        elec_ids = []
        for i, head in enumerate(nsx_file.extended_headers):
            elec_ids.append(head['ElectrodeID'])

    # Load data for all of the electrodes
    if verbose: print("Loading raw data from NSx file...", end='')
    chan_data = nsx_file.getdata(elec_ids=elec_ids)
    if verbose: print("done.")

    # Get resolution of signal in V
    if verbose: print("Processing electrode information.")
    cur_elec_ext_header_info = nsx_file.extended_headers[chan_data['ExtendedHeaderIndices'][0]]
    V_resolution = 0.001 * (
            cur_elec_ext_header_info['MaxAnalogValue'] - cur_elec_ext_header_info['MinAnalogValue']) / (
                           cur_elec_ext_header_info['MaxDigitalValue'] - cur_elec_ext_header_info[
                       'MinDigitalValue'])

    # Create electrode table for new electrodes
    elec_table_ind = np.ones((len(elec_ids), 1)) * np.NaN
    for cur_elec_ind, cur_elec in enumerate(elec_ids):
        # In order to create electrode table, we have to have indexes for each channel into electrode list in NWB file
        cur_elec_table_ind = np.where(nwb_file.electrodes[:, 0] == cur_elec)[0]
        if len(cur_elec_table_ind) == 0:
            raise ValueError("Couldn't find electrode %d in NWB file list." % (cur_elec))
        elec_table_ind[cur_elec_ind] = cur_elec_table_ind
    elec_table_ind = elec_table_ind.transpose().tolist()[0]  # Convert to list for create_electrode_table_region
    electrode_table_region = nwb_file.create_electrode_table_region(elec_table_ind, "Electrodes %s" % (elec_ids))

    # Write raw data for the electrode
    if verbose: print("\tAdding to NWB...", end='')
    cur_elec_ext_header_info = nsx_file.extended_headers[chan_data['ExtendedHeaderIndices'][0]]
    ephys_ts = ElectricalSeries('Raw Data for Channels %s' % (elec_ids),
                                chan_data['data'].transpose(),
                                electrode_table_region,
                                starting_time=chan_data['start_time_s'],
                                rate=chan_data['samp_per_s'],
                                resolution=V_resolution,
                                conversion=float(0.001),
                                # this is to convert the raw data (in mV) to V as expected by NWB
                                comments="Filter: High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)" % (
                                    cur_elec_ext_header_info['HighFreqCorner'],
                                    cur_elec_ext_header_info['HighFreqOrder'],
                                    cur_elec_ext_header_info['HighFreqType'],
                                    cur_elec_ext_header_info['LowFreqCorner'],
                                    cur_elec_ext_header_info['LowFreqOrder'],
                                    cur_elec_ext_header_info['LowFreqType']),
                                description="Raw data from NS6 file.")
    nwb_file.add_acquisition(ephys_ts)
    if verbose: print("done")

    # Write the file
    if verbose: print("\tWriting NWB file and closing.")
    nwb_io.write(nwb_file)
    nwb_io.close()

    # Close the nsx file
    nsx_file.datafile.close()

    return nwb_file_name


def AddBlackrockLFPDataToNWB(blk_file_name, nwb_file_name='', verbose=False, elec_ids=None):
    """
        Copies filtered LFP data from the specified Blackrock file to Neurodata Without Borders (NWB) file.
        User must specify the file with the LFP data (at desired sampling rate and with desired filtering).
        Multiple calls can be used to load multiple LFP data from multiple files.

        Typically, the NWB file will already be initialized by calling InitializeNWBFromBlackrock.
        Multiple electrodes can be added at once or with separate calls to this function.

        :param blk_file_name: {str} full path of Blackrock file to convert to NWB. If empty, will open dialog.
        :param nwb_file_name: [optional] {str} full path of NWB file to export to. Default is to change blackrock extension to 'nwb'
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :param elec_ids: [optional] {list} List of electrode IDs to copy over. If empty, all are copied
        :return: {str} filename of NWB file (empty if error)
        """
    # Check to see if user specified a Blackrock filename
    if not blk_file_name:  # no file name passed
        # Ask user to specify a file
        if 'app' not in locals():
            app = QApplication([])
        blk_file_name = QFileDialog.getOpenFileName(QFileDialog(), "Select File", getcwd(),
                                                    "Blackrock Data File (*.ns*)")
        blk_file_name = blk_file_name[0]

    # Check to see if valid nwb_file_name is passed
    if not nwb_file_name:
        nwb_file_name = path.splitext(blk_file_name)[0] + '.nwb'
    if verbose: print("Writing to NWB data file %s" % (nwb_file_name))

    # Initialize the Blackrock file
    try:
        # Requires the raw data to be imported
        nsx_file = brpylib.NsxFile(blk_file_name)
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        raise FileNotFoundError("Couldn't open Blackrock file. Error: %s" % e)

    # Initialize the NWB file
    nwb_io = []
    try:
        if not path.isfile(nwb_file_name):
            # Initialize NWB file
            if verbose: print("NWB file doesn't exist. Creating new one: %s..." % (nwb_file_name))
            InitializeNWBFromBlackrock(blk_file_name, nwb_file_name, verbose=True)

        # Append to existing file
        if verbose: print("Opening NWB file %s..." % (nwb_file_name), end='')
        nwb_file_append = True
        nwb_io = NWBHDF5IO(nwb_file_name, mode='a')
        nwb_file = nwb_io.read()
        if verbose: print("done.")
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        if nwb_io: nwb_io.close()
        nsx_file.datafile.close()
        raise FileExistsError("Couldn't open NWB file. Error: %s" % e)

    # Validate the elec_ids list
    if not elec_ids:
        # Grab all of the data
        elec_ids = []
        for i, head in enumerate(nsx_file.extended_headers):
            elec_ids.append(head['ElectrodeID'])

    # Load data for all of the electrodes
    if verbose: print("Loading LFP data from NSx file...", end='')
    chan_data = nsx_file.getdata(elec_ids=elec_ids)
    if verbose: print("done.")

    # Get resolution of signal in V
    if verbose: print("Processing electrode information.")
    cur_elec_ext_header_info = nsx_file.extended_headers[chan_data['ExtendedHeaderIndices'][0]]
    V_resolution = 0.001 * (
            cur_elec_ext_header_info['MaxAnalogValue'] - cur_elec_ext_header_info['MinAnalogValue']) / (
                           cur_elec_ext_header_info['MaxDigitalValue'] - cur_elec_ext_header_info[
                       'MinDigitalValue'])

    # Create electrode table for new electrodes
    elec_table_ind = np.ones((len(elec_ids), 1)) * np.NaN
    for cur_elec_ind, cur_elec in enumerate(elec_ids):
        # In order to create electrode table, we have to have indexes for each channel into electrode list in NWB file
        cur_elec_table_ind = np.where(nwb_file.electrodes[:, 0] == cur_elec)[0]
        if len(cur_elec_table_ind) == 0:
            raise ValueError("Couldn't find electrode %d in NWB file list." % (cur_elec))
        elec_table_ind[cur_elec_ind] = cur_elec_table_ind
    elec_table_ind = elec_table_ind.transpose().tolist()[0]  # Convert to list for create_electrode_table_region
    if verbose: print("Creating electrode table for electrodes %s" % (elec_ids))
    electrode_table_region = nwb_file.create_electrode_table_region(elec_table_ind, "Electrodes %s" % (elec_ids))

    # Create processing module for saving the data
    cur_elec_ext_header_info = nsx_file.extended_headers[chan_data['ExtendedHeaderIndices'][0]]
    if 'ecephys' not in nwb_file.processing.keys():
        if verbose: print("Processing module for extracellular electrophysiology (ecephys) does not exist. Creating.")
        ecephys_module = ProcessingModule(name='ecephys',
                                          description="Processing module for recordings from %s." % (
                                          path.split(blk_file_name)[1]))
        nwb_file.add_processing_module(ecephys_module)

    # Create data interface for the filtered LFP data
    cur_lfp = FilteredEphys(name="LFP data from %s" % (path.split(blk_file_name)[1]))
    if verbose: print("Created FilteredEphys data interface (%s)" % (cur_lfp.name))

    # Write LFP data for the electrode
    if verbose: print("\tAdding to NWB...", end='')
    cur_elec_ext_header_info = nsx_file.extended_headers[chan_data['ExtendedHeaderIndices'][0]]
    ephys_ts = ElectricalSeries('LFP Data for Channels %s' % (elec_ids),
                                chan_data['data'].transpose(),
                                electrode_table_region,
                                starting_time=chan_data['start_time_s'],
                                rate=chan_data['samp_per_s'],
                                resolution=V_resolution,
                                conversion=float(0.001),
                                # this is to convert the raw data (in mV) to V as expected by NWB
                                description="LFP data from %s. Filter: High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)" % (
                                    path.split(blk_file_name)[1],
                                    cur_elec_ext_header_info['HighFreqCorner'],
                                    cur_elec_ext_header_info['HighFreqOrder'],
                                    cur_elec_ext_header_info['HighFreqType'],
                                    cur_elec_ext_header_info['LowFreqCorner'],
                                    cur_elec_ext_header_info['LowFreqOrder'],
                                    cur_elec_ext_header_info['LowFreqType']))
    cur_lfp.add_electrical_series(ephys_ts)
    if verbose: print("done.")

    # Add LFP container to NWB file
    if verbose: print("\tAdding LFP data to processing module.")
    nwb_file.processing['ecephys'].add(cur_lfp)

    # Write the file
    if verbose: print("\tWriting NWB file and closing.")
    nwb_io.write(nwb_file)
    nwb_io.close()

    # Close the nsx file
    nsx_file.datafile.close()

    return nwb_file_name


def AddBlackrockSpikeDataToNWB(blk_nev_file_name, nwb_file_name='', nwb_file=None, elec_ids=None, verbose=False, add_units=False,
                               unit_info=None):
    """
        Copies filtered LFP data from the specified Blackrock file to Neurodata Without Borders (NWB) file.
        User must specify the file with the LFP data (at desired sampling rate and with desired filtering).
        Multiple calls can be used to load multiple LFP data from multiple files.

        Typically, the NWB file will already be initialized by calling InitializeNWBFromBlackrock.
        Multiple electrodes can be added at once or with separate calls to this function.

        :param blk_file_name: {str} full path of Blackrock file to convert to NWB. If empty, will open dialog.
        :param nwb_file_name: [optional] {str} full path of NWB file to export to. Default is to change blackrock extension to 'nwb'
        :param elec_ids: [optional] {list} List of electrode IDs to copy over. If empty, all are copied
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :param add_units: [optional] {bool} whether to add waveform labels as units in NWB file. Default is false.
        :return: {str} filename of NWB file (empty if error)
        """
    # Check to see if user specified a Blackrock filename
    if not blk_nev_file_name:  # no file name passed
        # Ask user to specify a file
        if 'app' not in locals():
            app = QApplication([])
        blk_nev_file_name = QFileDialog.getOpenFileName(QFileDialog(), "Select File", getcwd(),
                                                        "Blackrock Event Data File (*.nev)")
        blk_nev_file_name = blk_nev_file_name[0]

    # Initialize the Blackrock file
    try:
        # Requires the raw data to be imported
        nev_file = brpylib.NevFile(blk_nev_file_name)
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        raise FileNotFoundError("Couldn't open Blackrock file. Error: %s" % e)

    # Check to see if an open NWB file was passed, if so use that
    nwb_io = None
    if nwb_file == None:
        # Check to see if valid nwb_file_name is passed
        if not nwb_file_name:
            nwb_file_name = path.splitext(blk_nev_file_name)[0] + '.nwb'
        if verbose: print("Writing to NWB data file %s" % (nwb_file_name))

        # Initialize the NWB file
        try:
            if not path.isfile(nwb_file_name):
                # Initialize NWB file
                if verbose: print("NWB file doesn't exist. Creating new one: %s..." % (nwb_file_name))
                InitializeNWBFromBlackrock(blk_nev_file_name, nwb_file_name, verbose=True)

            # Append to existing file
            if verbose: print("Opening NWB file %s..." % (nwb_file_name), end='')
            nwb_file_append = True
            nwb_io = NWBHDF5IO(nwb_file_name, mode='a')
            nwb_file = nwb_io.read()
            if verbose: print("done.")
        except:  # catch *all* exceptions
            e = sys.exc_info()[0]
            if nwb_io: nwb_io.close()
            nev_file.datafile.close()
            raise FileExistsError("Couldn't open NWB file. Error: %s" % e)

    # Validate the elec_ids list
    if not elec_ids:
        # Grab all of the data
        elec_ids = []
        for i, head in enumerate(nev_file.extended_headers):
            if head['PacketID'] == 'NEUEVWAV':  # only take the waveform indexes for now
                elec_ids.append(head['ElectrodeID'])

    # Create processing module for saving ephys data
    if 'ecephys' not in nwb_file.processing.keys():
        if verbose: print("Processing module for extracellular electrophysiology (ecephys) does not exist. Creating.")
        ecephys_module = ProcessingModule(name='ecephys',
                                          description="Processing module for recordings from %s." % (
                                          path.split(blk_nev_file_name)[1]))
        nwb_file.add_processing_module(ecephys_module)

    # Create data interface for the spike waveforms
    cur_eventwaveform = EventWaveform(name="Spike waveforms from %s" % (path.split(blk_nev_file_name)[1]))
    if verbose:
        print("Created EventWaveform data interface (%s)" % (cur_eventwaveform.name))

    # Load data for all of the electrodes
    if verbose: print("Loading spike waveform data from NEV file...", end='')
    chan_data = nev_file.getdata(elec_ids=elec_ids)
    if verbose: print("done.")

    # Loop through and save spiking data for each electrode --  we are doing this separately for each electrode so that
    # we can save information about the threshold used for each channel in the description.
    if verbose: print("Writing spike waveforms data to NWB.")
    for cur_elec_ind, cur_elec in enumerate(elec_ids):
        # Find indices in the extended header information
        nev_wav_ind = chan_data['spike_events']['NEUEVWAV_HeaderIndices'][0]
        nev_flt_ind = 0
        while (nev_flt_ind < nev_file.basic_header['NumExtendedHeaders']) and \
                not ((cur_elec == nev_file.extended_headers[nev_flt_ind]['ElectrodeID']) and \
                     (nev_file.extended_headers[nev_flt_ind]['PacketID'] == 'NEUEVFLT')):
            nev_flt_ind = nev_flt_ind + 1
        # If we made it past the end, this electrode doesn't exist
        if (nev_wav_ind >= nev_file.basic_header['NumExtendedHeaders']) or \
                (nev_flt_ind >= nev_file.basic_header['NumExtendedHeaders']):
            print("Couldn't find specified electrode %d in NSx file. Skipping." % (cur_elec))
            continue
        if verbose: print("\tWriting spike waveforms for electrode %d." % (cur_elec))

        # Create electrode table
        cur_elec_table_ind = np.where(nwb_file.electrodes[:, 0] == cur_elec)[0]
        if len(cur_elec_table_ind) == 0:
            raise ValueError("Couldn't find electrode %d in NWB file list." % (cur_elec))
        electrode_table_region = nwb_file.create_electrode_table_region([cur_elec_table_ind],
                                                                        "Electrode %s" % (cur_elec))

        # Write raw data for the electrode
        if verbose: print("\tAdding to NWB...", end='')
        # Create timestamp array
        ts = np.array(chan_data['spike_events']['TimeStamps'][cur_elec_ind]) / nev_file.basic_header[
            'SampleTimeResolution']
        # Create spike event series
        spike_wf = SpikeEventSeries('Spike waveforms for Channel %s' % (cur_elec),
                                    chan_data['spike_events']['Waveforms'][cur_elec_ind],
                                    ts,
                                    electrode_table_region,
                                    conversion=1.0 * 10 ** -9,
                                    # these are stored in nV, so they need to be converted to V
                                    description="Spikes waveforms from channel %d in NEV file. Spikes were collected with threshold: Low: %d, Energy: %d, High: %d and filter: High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)" % (
                                        cur_elec,
                                        nev_file.extended_headers[nev_wav_ind]['LowThreshold'],
                                        nev_file.extended_headers[nev_wav_ind]['EnergyThreshold'],
                                        nev_file.extended_headers[nev_wav_ind]['HighThreshold'],
                                        nev_file.extended_headers[nev_flt_ind]['HighFreqCorner'],
                                        nev_file.extended_headers[nev_flt_ind]['HighFreqOrder'],
                                        nev_file.extended_headers[nev_flt_ind]['HighFreqType'],
                                        nev_file.extended_headers[nev_flt_ind]['LowFreqCorner'],
                                        nev_file.extended_headers[nev_flt_ind]['LowFreqOrder'],
                                        nev_file.extended_headers[nev_flt_ind]['LowFreqType']))

        # Add to our event waveform module
        cur_eventwaveform.add_spike_event_series(spike_wf)
        if verbose: print("done")

    # Add spike waveform container to NWB file
    if verbose: print("\tAdding spike waveform data to processing module.")
    nwb_file.processing['ecephys'].add(cur_eventwaveform)

    if add_units:
        # Check to see if user rating column already exists
        if (not nwb_file.units) or ('UserRating' not in nwb_file.units.colnames):
            nwb_file.add_unit_column('UserRating', 'Quality of isolation, as rated by user.')
        if (not nwb_file.units) or ('UserComment' not in nwb_file.units.colnames):
            nwb_file.add_unit_column('UserComment', 'Comments on spike isolation from user.')

        # Add units to the NWB file
        if verbose: print("Adding units to NWB file.")
        for cur_elec_ind, cur_elec in enumerate(elec_ids):
            # Determine the electrode index in NWB file for this electrode id
            cur_elec_table_ind = 0
            while (cur_elec_table_ind < len(nwb_file.electrodes)) and (
                    nwb_file.electrodes[cur_elec_table_ind, 0] != cur_elec):
                cur_elec_table_ind = cur_elec_table_ind + 1
            if cur_elec_table_ind >= len(nwb_file.electrodes):
                raise ValueError("Couldn't find electrode %d in NWB file list." % (cur_elec))

            # Get timestamps for this electrode
            elec_ts = np.array(chan_data['spike_events']['TimeStamps'][cur_elec_ind]) / nev_file.basic_header[
                'SampleTimeResolution']

            # Get list of unique waveform classifications
            uniq_class = np.unique(chan_data['spike_events']['Classification'][cur_elec_ind])
            if verbose: print("\tElectrode %d, found %d unique units (%s)." % (cur_elec, len(uniq_class), uniq_class))

            # Loop through unique classes and add units
            for cur_class in uniq_class:
                # Find the index in the unit_info list
                if not unit_info:
                    cur_unit_info = {'UserRating': -1, 'UserComment': ''}
                else:
                    i = 0
                    while (i < len(unit_info)) and (unit_info[i]['electrode'] != cur_elec) and (
                            unit_info[i]['unit'] != cur_class):
                        i = i + 1
                    if i >= len(unit_info):
                        cur_unit_info = {'UserRating': -1, 'UserComment': ''}
                    else:
                        cur_unit_info = unit_info[i]

                # Extract waveforms and timestamps
                class_ind = (np.array(chan_data['spike_events']['Classification'][cur_elec_ind]) == cur_class)
                wf = chan_data['spike_events']['Waveforms'][cur_elec_ind]
                ts = elec_ts[class_ind]

                # Calculate the observation interval
                obs_intervals = np.zeros((1, 2))
                obs_intervals[0, 0] = np.min(ts)
                obs_intervals[0, 1] = np.max(ts)

                # Add to NWB file
                nwb_file.add_unit(spike_times=ts,
                                  obs_intervals=obs_intervals,
                                  electrodes=[cur_elec_table_ind],
                                  waveform_mean=np.nanmean(wf, axis=0),
                                  waveform_sd=np.nanstd(wf, axis=0),
                                  UserRating=cur_unit_info['UserRating'],
                                  UserComment=cur_unit_info['UserComment'])
                if verbose: print("\t\tAdded class %s." % (cur_class))

    # Close the nev file
    nev_file.datafile.close()

    # Did we open a file?  Then close it...
    if nwb_io != None:
        # Write the file
        if verbose: print("\tWriting NWB file and closing.")
        nwb_io.write(nwb_file)
        nwb_io.close()
        return nwb_file_name
    else:
        return nwb_file


def AddBlackrockAnalogDataToNWB(blk_file_name, nwb_file_name='', signal_info=None, module_name='behavior',
                                verbose=False):
    """
        Copies analog (continuous) data from the specified Blackrock file to Neurodata Without Borders (NWB) file.
        This is usually continuous signals about behavior (joystick position, screen refreshes, etc).
        User should provide information about the signals to help end users understand what is in each signal.
        Multiple calls can be used to load multiple data signals from multiple files.

        Typically, the NWB file will already be initialized by calling InitializeNWBFromBlackrock.

        :param blk_file_name: {str} full path of Blackrock file to convert to NWB. If empty, will open dialog.
        :param nwb_file_name: [optional] {str} full path of NWB file to export to. Default is to change blackrock extension to 'nwb'
        :param signal_info: [optional] {list} List of dictionaries with information about the signals to save.
        :param module_name: [optional] {str} Name of module to store data.  Usually 'behavior' but could also be 'ecephys' or 'misc'
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :return: {str} filename of NWB file (empty if error)
        """
    # Check to see if user specified a Blackrock filename
    if not blk_file_name:  # no file name passed
        # Ask user to specify a file
        if 'app' not in locals():
            app = QApplication([])
        blk_file_name = QFileDialog.getOpenFileName(QFileDialog(), "Select File", getcwd(),
                                                    "Blackrock Data File (*.ns*)")
        blk_file_name = blk_file_name[0]

    # Check to see if valid nwb_file_name is passed
    if not nwb_file_name:
        nwb_file_name = path.splitext(blk_file_name)[0] + '.nwb'
    if verbose: print("Writing to NWB data file %s" % (nwb_file_name))

    # Initialize the Blackrock file
    try:
        nsx_file = brpylib.NsxFile(blk_file_name)
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        raise FileNotFoundError("Couldn't open Blackrock file. Error: %s" % e)

    # Initialize the NWB file
    nwb_io = []
    try:
        if not path.isfile(nwb_file_name):
            # Initialize NWB file
            if verbose: print("NWB file doesn't exist. Creating new one: %s..." % (nwb_file_name))
            InitializeNWBFromBlackrock(blk_file_name, nwb_file_name, verbose=True)

        # Append to existing file
        if verbose: print("Opening NWB file %s..." % (nwb_file_name), end='')
        nwb_file_append = True
        nwb_io = NWBHDF5IO(nwb_file_name, mode='a')
        nwb_file = nwb_io.read()
        if verbose: print("done.")
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        if nwb_io: nwb_io.close()
        nsx_file.datafile.close()
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
            nsx_ind = 0
            while (nsx_ind < len(nsx_file.extended_headers)) and (
                    nsx_file.extended_headers[nsx_ind]['ElectrodeLabel'] != cur_signal_info['label']):
                nsx_ind = nsx_ind + 1
            if nsx_ind >= len(nsx_file.extended_headers):
                raise ValueError("Couldn't find the specified signal %s" % (cur_signal_info['label']))
            signal_info[cur_signal_ind]['elec_id'] = nsx_file.extended_headers[nsx_ind]['ElectrodeID']
        elec_ids.append(signal_info[cur_signal_ind]['elec_id'])

    # Load data for all of the electrodes
    if verbose: print("Loading analog data from NSx file...", end='')
    chan_data = nsx_file.getdata(elec_ids=elec_ids)
    if verbose: print("done.")

    # Create processing module for saving data
    if module_name not in nwb_file.processing.keys():
        if verbose: print("Specified processing module (%s) does not exist. Creating." % (module_name))
        signal_module = ProcessingModule(name=module_name,
                                         description="Processing module for continuous signal data from %s." % (
                                             path.split(blk_file_name)[1]))
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
        # Get extended header information
        cur_elec_ext_header_info = nsx_file.extended_headers[chan_data['ExtendedHeaderIndices'][cur_signal_ind]]

        # Calculate resolution of signal in V
        V_resolution = 0.001 * (
                cur_elec_ext_header_info['MaxAnalogValue'] - cur_elec_ext_header_info['MinAnalogValue']) / (
                               cur_elec_ext_header_info['MaxDigitalValue'] - cur_elec_ext_header_info[
                           'MinDigitalValue'])

        # Create time series
        cur_data_interface.create_timeseries(name=cur_signal_info['name'],
                                             comments=cur_signal_info['comments'],
                                             unit="V",
                                             resolution=V_resolution,
                                             conversion=0.001,
                                             starting_time=chan_data['start_time_s'],
                                             rate=chan_data['samp_per_s'],
                                             description="Signal %s (Electrode ID: %d, Physical Connector: %d, Connector Pin: %d) from %s. Filter: High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)" % (
                                                 cur_signal_info['label'],
                                                 cur_elec_ext_header_info['ElectrodeID'],
                                                 cur_elec_ext_header_info['PhysicalConnector'],
                                                 cur_elec_ext_header_info['ConnectorPin'],
                                                 path.split(blk_file_name)[1],
                                                 cur_elec_ext_header_info['HighFreqCorner'],
                                                 cur_elec_ext_header_info['HighFreqOrder'],
                                                 cur_elec_ext_header_info['HighFreqType'],
                                                 cur_elec_ext_header_info['LowFreqCorner'],
                                                 cur_elec_ext_header_info['LowFreqOrder'],
                                                 cur_elec_ext_header_info['LowFreqType']))
        if verbose: print("\tAdded %s." % (cur_signal_info['label']))

    # Add data interface to module in NWB file
    if verbose: print("Adding data interface to module.")
    nwb_file.processing[module_name].add(cur_data_interface)

    # Write the file
    if verbose: print("Writing NWB file and closing.")
    nwb_io.write(nwb_file)
    nwb_io.close()

    # Close the nsx file
    nsx_file.datafile.close()

    return nwb_file_name


def AddBlackrockEyeDataToNWB(blk_file_name, nwb_file_name='', x_elec_id=None, y_elec_id=None, x_y_conversion=1.0,
                             x_y_reference_frame='Fixation', pupil_elec_id=None, pupil_conversion=1.0, verbose=False):
    """
        Copies analog (continuous) data from the specified Blackrock file to Neurodata Without Borders (NWB) file.
        This is only for continuous signals about eye position -- user provides the filename and electrodes carrying
        information about eye position and/or pupil diameter.

        Typically, the NWB file will already be initialized by calling InitializeNWBFromBlackrock.

        :param blk_file_name: {str} full path of Blackrock file to convert to NWB. If empty, will open dialog.
        :param nwb_file_name: [optional] {str} full path of NWB file to export to. Default is to change blackrock extension to 'nwb'
        :param x_elec_id: [optional] {int} Electrode ID for input carrying the horizontal (x) position of the eye
        :param y_elec_id: [optional] {int} Electrode ID for input carrying the vertical (y) position of the eye
        :param x_y_conversion: [optional] {float} Conversion factor to convert recorded signal to meaningful value (usually pixels)
        :param x_y_conversion: [optional] {str} Written description of what (0,0) is in eye data.  Default is 'Fixation'.
        :param pupil_elec_id: [optional] {int} Electrode ID for input carrying the pupil diameter of the eye
        :param pupil_conversion: [optional] {float} Conversion factor to convert recorded signal to meaningful value (usually pixels)
        :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
        :return: {str} filename of NWB file (empty if error)
        """

    # Check to see if user specified a Blackrock filename
    if not blk_file_name:  # no file name passed
        # Ask user to specify a file
        if 'app' not in locals():
            app = QApplication([])
        blk_file_name = QFileDialog.getOpenFileName(QFileDialog(), "Select File", getcwd(),
                                                    "Blackrock Data File (*.ns*)")
        blk_file_name = blk_file_name[0]

    # Check to see if valid nwb_file_name is passed
    if not nwb_file_name:
        nwb_file_name = path.splitext(blk_file_name)[0] + '.nwb'
    if verbose: print("Writing to NWB data file %s" % (nwb_file_name))

    # Initialize the Blackrock file
    try:
        nsx_file = brpylib.NsxFile(blk_file_name)
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        raise FileNotFoundError("Couldn't open Blackrock file. Error: %s" % e)

    # Initialize the NWB file
    nwb_io = []
    try:
        if not path.isfile(nwb_file_name):
            # Initialize NWB file
            if verbose: print("NWB file doesn't exist. Creating new one: %s..." % (nwb_file_name))
            InitializeNWBFromBlackrock(blk_file_name, nwb_file_name, verbose=verbose)

        # Append to existing file
        if verbose: print("Opening NWB file %s..." % (nwb_file_name), end='')
        nwb_file_append = True
        nwb_io = NWBHDF5IO(nwb_file_name, mode='a')
        nwb_file = nwb_io.read()
        if verbose: print("done.")
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        if nwb_io: nwb_io.close()
        nsx_file.datafile.close()
        raise FileExistsError("Couldn't open NWB file. Error: %s" % e)

    # Create processing module for saving data
    if 'behavior' not in nwb_file.processing.keys():
        if verbose: print("Behavior processing module does not exist. Creating.")
        signal_module = ProcessingModule(name='behavior',
                                         description="Processing module for behavior signals from %s." % (
                                             path.split(blk_file_name)[1]))
        nwb_file.add_processing_module(signal_module)

    # Add eye position data
    if x_elec_id or y_elec_id:
        if verbose: print("Adding eye position data...")

        # Check to see if labels were passed instead of electrode id #s
        if type(x_elec_id) is str:
            nsx_ind = 0
            while (nsx_ind < len(nsx_file.extended_headers)) and (
                    nsx_file.extended_headers[nsx_ind]['ElectrodeLabel'] != x_elec_id):
                nsx_ind = nsx_ind + 1
            if nsx_ind >= len(nsx_file.extended_headers):
                raise ValueError("Couldn't find the specified signal %s" % (x_elec_id))
            x_elec_id = nsx_file.extended_headers[nsx_ind]['ElectrodeID']
        if type(y_elec_id) is str:
            nsx_ind = 0
            while (nsx_ind < len(nsx_file.extended_headers)) and (
                    nsx_file.extended_headers[nsx_ind]['ElectrodeLabel'] != y_elec_id):
                nsx_ind = nsx_ind + 1
            if nsx_ind >= len(nsx_file.extended_headers):
                raise ValueError("Couldn't find the specified signal %s" % (y_elec_id))
            y_elec_id = nsx_file.extended_headers[nsx_ind]['ElectrodeID']
        if verbose: print("\tEye position data on electrodes %d and %d." % (x_elec_id, y_elec_id))

        # Create interface for eye data
        eye_data_interface = EyeTracking()

        # Load data for eye position
        if verbose: print("\tLoading eye data from NSx file...", end='')
        chan_data = nsx_file.getdata(elec_ids=[x_elec_id, y_elec_id])
        if verbose: print("done.")
        if not chan_data:
            raise ValueError("None of the electrode ids for eye position were found in the specified file.")

        # Create spatial series with eye data
        if verbose: print("\tCreating spatial series.")
        eye_data_interface.create_spatial_series("EyePosition from channels (%d, %d)" % (x_elec_id, y_elec_id),
                                                 chan_data['data'].transpose(),
                                                 conversion=x_y_conversion,
                                                 reference_frame=x_y_reference_frame,
                                                 starting_time=chan_data['start_time_s'],
                                                 rate=chan_data['samp_per_s'],
                                                 description='Eye tracking data recorded with Blackrock on channels %d and %d. Conversion factor should move it into degrees of visual angle.' % (
                                                 x_elec_id, y_elec_id))

        # Add data interface to module in NWB file
        if verbose: print("\tAdding eye position data interface to module.")
        nwb_file.processing['behavior'].add(eye_data_interface)

    # Add pupil diameter information
    if pupil_elec_id:
        if verbose: print("Adding pupil data...")

        # Check to see if labels were passed instead of electrode id #s
        if type(pupil_elec_id) is str:
            nsx_ind = 0
            while (nsx_ind < len(nsx_file.extended_headers)) and (nsx_file.extended_headers[nsx_ind]['ElectrodeLabel'] != pupil_elec_id):
                nsx_ind = nsx_ind + 1
            if nsx_ind >= len(nsx_file.extended_headers):
                raise ValueError("Couldn't find the specified signal %s" % (pupil_elec_id))
            pupil_elec_id = nsx_file.extended_headers[nsx_ind]['ElectrodeID']
        if verbose: print("\tPupil data on electrode %d." % (pupil_elec_id))

        # Create interface for eye data
        pupil_data_interface = PupilTracking()

        # Load data for pupil diameter
        if verbose: print("\tLoading pupil data from NSx file...", end='')
        chan_data = nsx_file.getdata(elec_ids=[pupil_elec_id])
        if verbose: print("done.")
        if not chan_data:
            raise ValueError("The electrode id for pupil diameter was not found in the specified file.")

        # Create spatial series with eye data
        if verbose: print("\tCreating time series.")
        pupil_data_interface.create_timeseries("EyePosition from channels (%d)" % (pupil_elec_id),
                                               chan_data['data'].transpose(),
                                               unit='pixel',
                                               conversion=pupil_conversion,
                                               starting_time=chan_data['start_time_s'],
                                               rate=chan_data['samp_per_s'],
                                               description='Pupil diamater data recorded with Blackrock on channel %d. Conversion factor should move it into pixels.' % (pupil_elec_id))

        # Add data interface to module in NWB file
        if verbose: print("\tAdding pupil diameter data interface to module.")
        nwb_file.processing['behavior'].add(pupil_data_interface)

    # Write the file
    if verbose: print("Writing NWB file and closing.")
    nwb_io.write(nwb_file)
    nwb_io.close()

    # Close the nsx file
    nsx_file.datafile.close()

    return nwb_file_name


def ConvertBlackrockEphysToNWB(blk_file_name, nwb_file_name='', verbose=False, experiment=None, electrode_info=None,
                               electrode_group=None, overwrite=False):
    """
    Copies electrophysiological data from Blackrock format to Neurodata Without Borders (NWB) format.
    User can provide information about experiment and electrodes to add to the NWB file.  If file already
    exists, data will be appended.

    :param blk_file_name: {str} full path of Blackrock file to convert to NWB. If empty, will open dialog.
    :param nwb_file_name: [optional] {str} full path of NWB file to export to. Default is to change blackrock extension to 'nwb'
    :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
    :param experiment: [optional] {dict} Dictionary that can contain information about experimenter, name, animal, lab, and institution
    :param electrode_info: [optional] {list} List of dictionaries carrying information about electrodes
    :param electrode_group: [optional] {list} List of dictionaries carrying information about electrode groupings
    :param overwrite: [optional] {bool} Whether to overwrite the NWB file if it already exists
    :return: {str} filename of NWB file (empty if error)
    """
    # Check to see if user specified a Blackrock filename
    if not blk_file_name:  # no file name passed
        # Ask user to specify a file
        if 'app' not in locals():
            app = QApplication([])
        blk_file_name = QFileDialog.getOpenFileName(QFileDialog(), "Select File", getcwd(),
                                                    "Blackrock Data File (*.ns6)")
        blk_file_name = blk_file_name[0]

    # Check to see if Blackrock files exist
    # NSx: Neural Stream X, where X indicates the sampling frequency of the continuous file. This data type contains
    #   continuous streamed data from the NSP.Depending on the sampling frequency, a second of data may contain up
    #   to 30, 000 samples in this file.
    #   NS1: Data sampled at 500 Hz
    blk_ns1_file_name = path.splitext(blk_file_name)[0] + '.ns1'
    if not path.isfile(blk_ns1_file_name):
        blk_ns1_file_name = ''
    #   NS2: Data sampled at 1 kHz
    blk_ns2_file_name = path.splitext(blk_file_name)[0] + '.ns2'
    if not path.isfile(blk_ns2_file_name):
        blk_ns2_file_name = ''
    #   NS3: Data sampled at 2 kHz
    blk_ns3_file_name = path.splitext(blk_file_name)[0] + '.ns3'
    if not path.isfile(blk_ns3_file_name):
        blk_ns3_file_name = ''
    #   NS4: Data sampled at 10 kHz
    blk_ns4_file_name = path.splitext(blk_file_name)[0] + '.ns4'
    if not path.isfile(blk_ns4_file_name):
        blk_ns4_file_name = ''
    #   NS5: Data sampled at 30 kHz
    blk_ns5_file_name = path.splitext(blk_file_name)[0] + '.ns5'
    if not path.isfile(blk_ns5_file_name):
        blk_ns5_file_name = ''

    # Check to see if valid nwb_file_name is passed
    if not nwb_file_name:
        nwb_file_name = path.splitext(blk_file_name)[0] + '.nwb'
    if verbose: print("Writing to NWB data file %s" % (nwb_file_name))

    # Initialize the Blackrock file
    try:
        # Requires the raw data to be imported
        nsx_file = brpylib.NsxFile(blk_file_name)
        # Process all of the NS* files
        nsn_file = []
        nsn_file_names = []
        if blk_ns1_file_name:
            nsn_file.append(brpylib.NsxFile(blk_ns1_file_name))
            nsn_file_names.append(blk_ns1_file_name)
        if blk_ns2_file_name:
            nsn_file.append(brpylib.NsxFile(blk_ns2_file_name))
            nsn_file_names.append(blk_ns2_file_name)
        if blk_ns3_file_name:
            nsn_file.append(brpylib.NsxFile(blk_ns3_file_name))
            nsn_file_names.append(blk_ns3_file_name)
        if blk_ns4_file_name:
            nsn_file.append(brpylib.NsxFile(blk_ns4_file_name))
            nsn_file_names.append(blk_ns4_file_name)
        if blk_ns5_file_name:
            nsn_file.append(brpylib.NsxFile(blk_ns5_file_name))
            nsn_file_names.append(blk_ns5_file_name)
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        raise FileNotFoundError("Couldn't open Blackrock file. Error: %s" % e)

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
                                  ', Date: ' + nsx_file.basic_header['TimeOrigin'].strftime(
            '%m/%d/%Y')  # description of the recording session
        if verbose: print("Experiment description: %s" % (nwb_session_description))
    # Define the NWB file identifier as the filename, as this should be unique
    nwb_identifier = path.split(blk_file_name)[1]  # unique file identifier

    # Initialize the NWB file
    try:
        if not overwrite and path.isfile(nwb_file_name):
            # Append to existing file
            if verbose: print("Opening existing NWB file %s..." % (nwb_file_name), end='')
            nwb_file_append = True
            nwb_io = NWBHDF5IO(nwb_file_name, mode='a')
            nwb_file = nwb_io.read()
            if verbose: print("done.")
        else:
            if verbose: print("Creating new NWB file %s..." % (nwb_file_name), end='')
            nwb_file_append = False
            nwb_io = NWBHDF5IO(nwb_file_name, mode='w')
            nwb_file = NWBFile(nwb_session_description, nwb_identifier, nsx_file.basic_header['TimeOrigin'],
                               file_create_date=datetime.now(tzlocal()),
                               experimenter=experimenter,
                               lab=lab,
                               institution=institution,
                               source_script='ConvertBlackrockToNWB.py',
                               source_script_file_name='ConvertBlackrockToNWB.py',
                               notes=nsx_file.basic_header['Comment'])
            if verbose: print("done.")
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        nwb_io.close()
        nsx_file.datafile.close()
        for cur_nsn_file_ind, cur_nsn_file in enumerate(nsn_file):
            cur_nsn_file.datafile.close()
        raise FileExistsError("Couldn't open NWB file. Error: %s" % e)
    if verbose: print("Created NWB file.")

    # Create device in the NWB file
    device = nwb_file.create_device(
        name='Blackrock (' + nsx_file.basic_header['FileTypeID'] + ', ' + nsx_file.basic_header['FileSpec'] + ')')
    if verbose: print("Created device.")

    # Make sure electrode input is provided and valid
    if (not electrode_info) or not (type(electrode_info) is list):
        if verbose: print("Creating electrode information from file.")
        # Create electrode information for each electrode
        electrode_info = []
        for cur_elec in range(nsx_file.basic_header['ChannelCount']):
            electrode_info.append({'x': 0.0, 'y': 0.0, 'z': 0.0, 'impedance': -1.0,
                                   'location': 'unknown', 'group': 'electrodes',
                                   'id': nsx_file.extended_headers[cur_elec]['ElectrodeID'],
                                   'filtering': 'High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)' %
                                                (nsx_file.extended_headers[0]['HighFreqCorner'],
                                                 nsx_file.extended_headers[0]['HighFreqOrder'],
                                                 nsx_file.extended_headers[0]['HighFreqType'],
                                                 nsx_file.extended_headers[0]['LowFreqCorner'],
                                                 nsx_file.extended_headers[0]['LowFreqOrder'],
                                                 nsx_file.extended_headers[0]['LowFreqType'])})
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
            raise ValueError(
                "Couldn't process electrode inputs. Must provide IDs for all electrodes or pass empty array.")
        # If filtering is not provided, add it from NSX file
        if ('filtering' not in cur_elec_info.keys()) or (not cur_elec_info['filtering']):
            electrode_info[cur_elec_info_ind]['filtering'] = 'High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)' % (
            nsx_file.extended_headers[0]['HighFreqCorner'], nsx_file.extended_headers[0]['HighFreqOrder'],
            nsx_file.extended_headers[0]['HighFreqType'], nsx_file.extended_headers[0]['LowFreqCorner'],
            nsx_file.extended_headers[0]['LowFreqOrder'], nsx_file.extended_headers[0]['LowFreqType'])
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
        if grp_elec_ind >= len(electrode_group): grp_elec_ind = -1
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

    # Loop through and save raw data for each electrode
    if verbose: print("Writing raw data to NWB.")
    for cur_elec_ind, cur_elec in enumerate(electrode_info):
        # Find matching electrode in nsx file
        nsx_elec_ind = 0
        while (nsx_elec_ind < nsx_file.basic_header['ChannelCount']) and (
                cur_elec['id'] != nsx_file.extended_headers[nsx_elec_ind]['ElectrodeID']):
            nsx_elec_ind = nsx_elec_ind + 1
        # If we made it past the end, this electrode doesn't exist
        if nsx_elec_ind >= nsx_file.basic_header['ChannelCount']:
            warnings.warn("Couldn't find specified electrode %d in NS6 file. Skipping." % (cur_elec['id']), UserWarning)
            continue
        if verbose: print("Writing raw data for electrode %d." % (cur_elec['id']))

        # Get raw data from NSX file
        if verbose: print("\tLoading data from NSx file.")
        chan_data = nsx_file.getdata(elec_ids=cur_elec['id'])

        # Get resolution of signal in V
        cur_elec_ext_header_info = nsx_file.extended_headers[chan_data['ExtendedHeaderIndices'][0]]
        V_resolution = 1000 * (
                    cur_elec_ext_header_info['MaxAnalogValue'] - cur_elec_ext_header_info['MinAnalogValue']) / (
                               cur_elec_ext_header_info['MaxDigitalValue'] - cur_elec_ext_header_info[
                           'MinDigitalValue'])

        # Create electrode table
        elec_table_ind = 0
        while (elec_table_ind < len(nwb_file_elec_list)) and (nwb_file_elec_list[elec_table_ind] != cur_elec['id']):
            elec_table_ind = elec_table_ind + 1
        if elec_table_ind >= len(nwb_file_elec_list):
            raise AssertionError(
                "Couldn't find electrode %d in NWB file list. This shouldn't happen." % (cur_elec['id']))
        electrode_table_region = nwb_file.create_electrode_table_region([elec_table_ind],
                                                                        "Electrode %d" % (elec_table_ind))

        # Write raw data for the electrode
        if verbose: print("\tAdding to NWB...", end='')
        ephys_ts = ElectricalSeries('Raw Data for Chan %d (%s, Physical Connector %d, Connector Pin: %d).' % (
            cur_elec['id'], cur_elec_ext_header_info['ElectrodeLabel'], cur_elec_ext_header_info['PhysicalConnector'],
            cur_elec_ext_header_info['ConnectorPin']),
                                    chan_data['data'][0],
                                    electrode_table_region,
                                    starting_time=chan_data['start_time_s'],
                                    rate=chan_data['samp_per_s'],
                                    resolution=V_resolution,
                                    conversion=float(1000),
                                    # this is to convert the raw data (in uV) to V as expected by NWB
                                    comments="Filter: High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)" % (
                                        cur_elec_ext_header_info['HighFreqCorner'],
                                        cur_elec_ext_header_info['HighFreqOrder'],
                                        cur_elec_ext_header_info['HighFreqType'],
                                        cur_elec_ext_header_info['LowFreqCorner'],
                                        cur_elec_ext_header_info['LowFreqOrder'],
                                        cur_elec_ext_header_info['LowFreqType']),
                                    description="Raw data from NS6 file.")
        nwb_file.add_acquisition(ephys_ts)
        if verbose: print("done")

    # Loop through and save LFP data (from NS1 through NS5) for each electrode
    for cur_nsn_file_ind, cur_nsn_file in enumerate(nsn_file):
        if verbose: print("Processing electrodes from %s." % (nsn_file_names[cur_nsn_file_ind]))
        cur_lfp = []
        for cur_elec_ind, cur_elec in enumerate(electrode_info):
            # Find matching electrode in nsx file
            nsn_elec_ind = 0
            while (nsn_elec_ind < len(cur_nsn_file.extended_headers)) and (
                    cur_elec['id'] != cur_nsn_file.extended_headers[nsn_elec_ind]['ElectrodeID']):
                nsn_elec_ind = nsn_elec_ind + 1
            # If we made it past the end, this electrode doesn't exist
            if nsn_elec_ind >= cur_nsn_file.basic_header['ChannelCount']:
                print("Couldn't find specified electrode %d in NSx file. Skipping." % (cur_elec['id']))
                continue
            if verbose: print("\tWriting LFP data for electrode %d." % (cur_elec['id']))

            # If we haven't yet created an LFP container in the NWB file, we should do so now
            if not cur_lfp:
                ecephys_module = ProcessingModule(name="LFP data from %s" % (nsn_file_names[cur_nsn_file_ind]),
                                                  description="Filter: High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)" % (
                                                      cur_elec_ext_header_info['HighFreqCorner'],
                                                      cur_elec_ext_header_info['HighFreqOrder'],
                                                      cur_elec_ext_header_info['HighFreqType'],
                                                      cur_elec_ext_header_info['LowFreqCorner'],
                                                      cur_elec_ext_header_info['LowFreqOrder'],
                                                      cur_elec_ext_header_info['LowFreqType']))
                nwb_file.add_processing_module(ecephys_module)
                cur_lfp = FilteredEphys(name="LFP data from %s" % (nsn_file_names[cur_nsn_file_ind]))

            # Get raw data from NSn file
            if verbose: print("\t\tLoading data from NSn file.")
            chan_data = cur_nsn_file.getdata(elec_ids=cur_elec['id'])

            # Get resolution of signal in V
            cur_elec_ext_header_info = cur_nsn_file.extended_headers[chan_data['ExtendedHeaderIndices'][0]]
            V_resolution = 1000 * (
                        cur_elec_ext_header_info['MaxAnalogValue'] - cur_elec_ext_header_info['MinAnalogValue']) / (
                                   cur_elec_ext_header_info['MaxDigitalValue'] - cur_elec_ext_header_info[
                               'MinDigitalValue'])

            # Create electrode table
            elec_table_ind = 0
            while (elec_table_ind < len(nwb_file_elec_list)) and (nwb_file_elec_list[elec_table_ind] != cur_elec['id']):
                elec_table_ind = elec_table_ind + 1
            if elec_table_ind >= len(nwb_file_elec_list):
                raise AssertionError(
                    "Couldn't find electrode %d in NWB file list. This shouldn't happen." % (cur_elec['id']))
            electrode_table_region = nwb_file.create_electrode_table_region([elec_table_ind],
                                                                            "Electrode %d" % (elec_table_ind))

            # Write raw data for the electrode
            if verbose: print("\tAdding to LFP...", end='')
            ephys_ts = ElectricalSeries('LFP Data for Chan %d (%s, Physical Connector %d, Connector Pin: %d).' % (
                cur_elec['id'], cur_elec_ext_header_info['ElectrodeLabel'],
                cur_elec_ext_header_info['PhysicalConnector'], cur_elec_ext_header_info['ConnectorPin']),
                                        chan_data['data'][0],
                                        electrode_table_region,
                                        starting_time=chan_data['start_time_s'],
                                        rate=chan_data['samp_per_s'],
                                        resolution=V_resolution,
                                        conversion=float(1000),
                                        # this is to convert the raw data (in uV) to V as expected by NWB
                                        comments="Filter: High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)" % (
                                            cur_elec_ext_header_info['HighFreqCorner'],
                                            cur_elec_ext_header_info['HighFreqOrder'],
                                            cur_elec_ext_header_info['HighFreqType'],
                                            cur_elec_ext_header_info['LowFreqCorner'],
                                            cur_elec_ext_header_info['LowFreqOrder'],
                                            cur_elec_ext_header_info['LowFreqType']),
                                        description="Raw data from NS6 file.")
            cur_lfp.add_electrical_series(ephys_ts)
            if verbose: print("done.")

        # Add LFP container to NWB file
        if cur_lfp:
            if verbose: print("\tAdding LFP data to module.")
            nwb_file.processing["LFP data from %s" % (nsn_file_names[cur_nsn_file_ind])].add(cur_lfp)

        # Close the NSn file
        if verbose: print("\tClosing NSn file.")
        cur_nsn_file.datafile.close()

    # Close the NSx file
    if verbose: print("\tClosing NSx file.")
    nsx_file.datafile.close()

    # Write the file
    if verbose: print("\tWriting NWB file and closing.")
    nwb_io.write(nwb_file)
    nwb_io.close()

    # Return the file name
    return nwb_file_name


def ConvertBlackrockSpikesToNWB(blk_file_name='', nwb_file_name='', verbose=False, experiment=None, electrode_info=None,
                                electrode_group=None, overwrite=False):
    """
    Copy spikes (and waveforms) from Blackrock data file to Neurodata Without Borders (NWB) format.
    If file already exists, can append NWB file with spikes.

    :param blk_file_name: {str} full path of Blackrock file to convert to NWB. If empty, will open dialog.
    :param nwb_file_name: [optional] {str} full path of NWB file to export to. Default is to change blackrock extension to 'nwb'
    :param verbose: [optional] {bool} whether to print updates while converting. Default is false.
    :param experiment: [optional] {dict} Dictionary that can contain information about experimenter, name, animal, lab, and institution
    :param electrode_info: [optional] {list} List of dictionaries carrying information about electrodes
    :param electrode_group: [optional] {list} List of dictionaries carrying information about electrode groupings
    :param overwrite: [optional] {bool} Whether to overwrite the NWB file if it already exists
    :return: {str} filename of NWB file (empty if error)
    """
    # Check to see if user specified a Blackrock filename
    if not blk_file_name:  # no file name passed
        # Ask user to specify a file
        if 'app' not in locals():
            app = QApplication([])
        blk_file_name = QFileDialog.getOpenFileName(QFileDialog(), "Select File", getcwd(),
                                                    "Blackrock Data File (*.ns6)")
        blk_file_name = blk_file_name[0]

    # Check to see if Blackrock files exist
    # NEV -- contains events recorded during an experimental session. These events may include information
    #   on threshold crossed spike waveforms, information received through the digital input and the serial input,
    #   tracking information recorded in NeuroMotive Video Tracking System, text comments, or other custom comments
    #   sent to the NSP during an experiment.
    blk_nev_file_name = path.splitext(blk_file_name)[0] + '.nev'
    if not path.isfile(blk_nev_file_name):
        if verbose: print("Couldn't find associated NEV file. Skipping.")
        blk_nev_file_name = ''

    # Check to see if valid nwb_file_name is passed
    if not nwb_file_name:
        nwb_file_name = path.splitext(blk_file_name)[0] + '.nwb'
    if verbose: print("Writing to NWB data file %s" % (nwb_file_name))

    # Initialize the Blackrock file
    try:
        # Requires the raw data to be imported
        nsx_file = brpylib.NsxFile(blk_file_name)
        if blk_nev_file_name:
            nev_file = brpylib.NevFile(blk_nev_file_name)
        else:
            raise FileNotFoundError("Couldn't find NEV file matching the specified Blackrock file.")
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        print("Couldn't open Blackrock file. Error: %s" % e)
        return False

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
                                  ', Date: ' + nsx_file.basic_header['TimeOrigin'].strftime(
            '%m/%d/%Y')  # description of the recording session
        if verbose: print("Experiment description: %s" % (nwb_session_description))
    # Define the NWB file identifier as the filename, as this should be unique
    nwb_identifier = path.split(blk_file_name)[1]  # unique file identifier

    # Initialize the NWB file
    try:
        if not overwrite and path.isfile(nwb_file_name):
            # Append to existing file
            if verbose: print("Opening existing NWB file %s..." % (nwb_file_name), end='')
            nwb_file_append = True
            nwb_io = NWBHDF5IO(nwb_file_name, mode='a')
            nwb_file = nwb_io.read()
            if verbose: print("done.")
        else:
            if verbose: print("Creating new NWB file %s..." % (nwb_file_name), end='')
            nwb_file_append = False
            nwb_io = NWBHDF5IO(nwb_file_name, mode='w')
            nwb_file = NWBFile(nwb_session_description, nwb_identifier, nsx_file.basic_header['TimeOrigin'],
                               file_create_date=datetime.now(tzlocal()),
                               experimenter=experimenter, lab=lab,
                               institution=institution,
                               source_script='ConvertBlackrockToNWB.py',
                               notes=nsx_file.basic_header['Comment'])
            if verbose: print("done.")
    except:  # catch *all* exceptions
        e = sys.exc_info()[0]
        raise FileExistsError("Couldn't open NWB file. Error: %s" % e)
    if verbose: print("Created NWB file.")

    # Create device in the NWB file
    device = nwb_file.create_device(
        name='Blackrock (' + nsx_file.basic_header['FileTypeID'] + ', ' + nsx_file.basic_header['FileSpec'] + ')')
    if verbose: print("Created device.")

    # Make sure electrode input is provided and valid
    if (not electrode_info) or not (type(electrode_info) is dict):
        if verbose: print("Creating electrode information from file.")
        # Create electrode information for each electrode
        electrode_info = []
        for cur_elec in range(nsx_file.basic_header['ChannelCount']):
            electrode_info.append({'x': 0.0, 'y': 0.0, 'z': 0.0, 'impedance': -1,
                                   'location': 'unknown', 'group': 'electrodes',
                                   'id': nsx_file.extended_headers[cur_elec]['ElectrodeID'],
                                   'filtering': 'High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)' %
                                                (nsx_file.extended_headers[0]['HighFreqCorner'],
                                                 nsx_file.extended_headers[0]['HighFreqOrder'],
                                                 nsx_file.extended_headers[0]['HighFreqType'],
                                                 nsx_file.extended_headers[0]['LowFreqCorner'],
                                                 nsx_file.extended_headers[0]['LowFreqOrder'],
                                                 nsx_file.extended_headers[0]['LowFreqType'])})
        if verbose: print("\tCreated %d electrodes." % (len(electrode_info)))

    # Fill in any missing information about electrodes
    ovr_elec_group_list = []  # will be used to keep track of electrode group list
    if verbose: print("Processing electrode information.")
    for cur_elec_info_ind, cur_elec_info in enumerate(electrode_info):
        # If id is not provided, then throw an error
        if ('id' not in cur_elec_info.keys()) or (not cur_elec_info['id']):
            raise ValueError(
                "Couldn't process electrode inputs. Must provide IDs for all electrodes or pass empty array.")
        # If filtering is not provided, add it from NSX file
        if ('filtering' not in cur_elec_info.keys()) or (not cur_elec_info['filtering']):
            cur_elec_info['filtering'] = 'High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)' % (
            nsx_file.extended_headers[0]['HighFreqCorner'], nsx_file.extended_headers[0]['HighFreqOrder'],
            nsx_file.extended_headers[0]['HighFreqType'], nsx_file.extended_headers[0]['LowFreqCorner'],
            nsx_file.extended_headers[0]['LowFreqOrder'], nsx_file.extended_headers[0]['LowFreqType'])
        # If other variables are not provided, just initialize them to defaults
        if ('x' not in cur_elec_info.keys()) or (not cur_elec_info['x']):
            cur_elec_info['x'] = 0.0
        if ('y' not in cur_elec_info.keys()) or (not cur_elec_info['y']):
            cur_elec_info['y'] = 0.0
        if ('z' not in cur_elec_info.keys()) or (not cur_elec_info['z']):
            cur_elec_info['z'] = 0.0
        if ('location' not in cur_elec_info.keys()) or (not cur_elec_info['location']):
            cur_elec_info['location'] = 'unknown'
        if ('impedance' not in cur_elec_info.keys()) or (not cur_elec_info['impedance']):
            cur_elec_info['impedance'] = -1
        if ('group' not in cur_elec_info.keys()) or (not cur_elec_info['group']):
            cur_elec_info['group'] = 'default electrode group'
        # Find the index to the electrode group
        grp_elec_ind = 0
        while (grp_elec_ind < len(electrode_group)) and (
                cur_elec_info['group'] != electrode_group[grp_elec_ind]['name']):
            grp_elec_ind = grp_elec_ind + 1
        # If we made it past the end, this electrode group doesn't exist, put it in the default
        if grp_elec_ind >= len(electrode_group): grp_elec_ind = -1
        # Save the index into group dictionary for easy reference
        cur_elec_info['electrode_group_index'] = grp_elec_ind
        ovr_elec_group_list.append(grp_elec_ind)

    # Create electrode groups
    nwb_electrode_group = []
    if verbose: print("Creating electrode groups.")
    for cur_elec_group in range(len(electrode_group) + 1):
        if (cur_elec_group - 1) in ovr_elec_group_list:
            if cur_elec_group == 0:
                # This is our default group for electrodes that don't fall into other groups
                nwb_electrode_group.append(nwb_file.create_electrode_group('default electrode group',
                                                                           description='Generic electrode group for ungrouped electrodes.',
                                                                           location='Unknown', device=device))
            else:
                # Create and add to our list
                nwb_electrode_group.append(nwb_file.create_electrode_group(electrode_group[cur_elec_group - 1]['name'],
                                                                           description=
                                                                           electrode_group[cur_elec_group - 1][
                                                                               'description'],
                                                                           location=electrode_group[cur_elec_group - 1][
                                                                               'location'], device=device))
        else:
            nwb_electrode_group.append([])
    if verbose: print("\tCreated %d electrode groups." % (len(electrode_group)))

    # Create electrodes in NWB file
    if verbose: print("Adding electrodes to NWB.")
    nwb_file_elec_list = []
    for cur_elec_ind, cur_elec in enumerate(electrode_info):
        # Add electrode to NWB file
        nwb_file.add_electrode(id=cur_elec['id'], x=cur_elec['x'], y=cur_elec['y'], z=cur_elec['z'],
                               imp=cur_elec['impedance'], location=cur_elec['location'],
                               filtering=cur_elec['filtering'],
                               group=nwb_electrode_group[cur_elec['electrode_group_index'] + 1])
        # Keep track of electrodes entered so we can index for the electrode table
        nwb_file_elec_list.append(cur_elec['id'])
    if verbose: print("\tAdded %d electrodes." % (len(electrode_info)))

    # Create EventWaveform module to hold data
    cur_eventwaveform = EventWaveform(name="Spike waveforms from %s." % (blk_file_name))

    # Loop through and save spiking data for each electrode
    if verbose: print("Writing raw data to NWB.")
    for cur_elec_ind, cur_elec in enumerate(electrode_info):
        # Find matching electrode in nsx file
        nev_wav_ind = 0
        while (cur_elec['id'] != nev_file.extended_headers[nev_wav_ind]['ElectrodeID']) and \
                (nev_wav_ind < nev_file.basic_header['NumExtendedHeaders']) and \
                (nev_file.extended_headers['PacketID'] == 'NEUEVWAV'):
            nev_wav_ind = nev_wav_ind + 1
        nev_flt_ind = 0
        while (cur_elec['id'] != nev_file.extended_headers[nev_flt_ind]['ElectrodeID']) and \
                (nev_flt_ind < nev_file.basic_header['NumExtendedHeaders']) and \
                (nev_file.extended_headers['PacketID'] == 'NEUEVFLT'):
            nev_flt_ind = nev_flt_ind + 1
        # If we made it past the end, this electrode doesn't exist
        if (nev_wav_ind >= nev_file.basic_header['NumExtendedHeaders']) or \
                (nev_flt_ind >= nev_file.basic_header['NumExtendedHeaders']):
            print("Couldn't find specified electrode %d in NSx file. Skipping." % (cur_elec['id']))
            continue
        if verbose: print("Writing spike waveforms for electrode %d." % (cur_elec['id']))

        # Get raw data from NSX file
        if verbose: print("\tLoading data from NEV file.")
        chan_data = nev_file.getdata(elec_ids=cur_elec['id'])
        if not chan_data:
            warnings.warn("Couldn't find specified electrode %d in NEV file. Skipping." % (cur_elec['id']), UserWarning)
            continue

        # Create electrode table
        elec_table_ind = 0
        while (elec_table_ind < len(nwb_file_elec_list)) and (nwb_file_elec_list[elec_table_ind] != cur_elec['id']):
            elec_table_ind = elec_table_ind + 1
        if elec_table_ind >= len(nwb_file_elec_list):
            raise AssertionError(
                "Couldn't find electrode %d in NWB file list. This shouldn't happen." % (cur_elec['id']))
        electrode_table_region = nwb_file.create_electrode_table_region(elec_table_ind,
                                                                        "Electrode %d" % (elec_table_ind))

        # Write raw data for the electrode
        if verbose: print("\tAdding to NWB...", end='')
        ts = np.array(chan_data['spike_events']['TimeStamps'][0]) / nev_file.basic_header['SampleTimeResolution']
        spike_wf = SpikeEventSeries('Spike waveforms for Chan %d (Physical Connector %d, Connector Pin: %d).' % (
            cur_elec['ElectrodeID'], cur_elec['PhysicalConnector'], cur_elec['ConnectorPin']),
                                    chan_data['spike_events']['Waveforms'][0],
                                    ts,
                                    electrode_table_region,
                                    conversion=1.0 * 10 ** -9,
                                    # these are stored in nV, so they need to be converted to V
                                    comments="Threshold: Low: %d, Energy: %d, High: %d; Filter: High Pass: %s (%d pole, %s), Low Pass: %s (%d pole, %s)" % (
                                        nev_file.extended_headers[nev_wav_ind]['LowThreshold'],
                                        nev_file.extended_headers[nev_wav_ind]['EnergyThreshold'],
                                        nev_file.extended_headers[nev_wav_ind]['HighThreshold'],
                                        nev_file.extended_headers[nev_flt_ind]['HighFreqCorner'],
                                        nev_file.extended_headers[nev_flt_ind]['HighFreqOrder'],
                                        nev_file.extended_headers[nev_flt_ind]['HighFreqType'],
                                        nev_file.extended_headers[nev_flt_ind]['LowFreqCorner'],
                                        nev_file.extended_headers[nev_flt_ind]['LowFreqOrder'],
                                        nev_file.extended_headers[nev_flt_ind]['LowFreqType']),
                                    description='Spike waveforms from Chan %d.' % (cur_elec['id']))

        # Add to our event waveform module
        cur_eventwaveform.add_spike_event_series(spike_wf)
        if verbose: print("done")

    # Close the NSx file
    nsx_file.datafile.close()

    # Write the file
    nwb_io.write(nwb_file)
    nwb_io.close()

    # Return the file name
    return nwb_file_name
