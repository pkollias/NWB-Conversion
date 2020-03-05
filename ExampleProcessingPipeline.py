import ConvertBlackrockToNWB
import ConvertPlexonToNWB
import NWBHelperFunctions
import NWBAnalysisFunctions
# import sys
import pandas
from datetime import datetime
import warnings

from pynwb import ProcessingModule
from pynwb.behavior import BehavioralEpochs
from pynwb.misc import IntervalSeries

from Gating.rec import *


def main():

    # init
    # args = sys.argv
    fsroot = '/Volumes'  # args[1]
    sess_iloc = 27  # int(args[2])
    md = MetaData(fsroot)
    db = md.db_base_loader(['sessions', 'trials', 'events', 'units', 'conditions', 'activity'])
    sessions, trials, events, units, conditions, activity =\
        db['sessions'], db['trials'], db['events'], db['units'], db['conditions'], db['activity']

    subject_to_info = {'Oscar': {'subject_id': 'Oscar', 'age': '8', 'description': 'NA',
                                 'genotype': 'NA', 'sex': 'M', 'species': 'Macaca mulatta', 'weight': '9.5',
                                 'date_of_birth': datetime.strptime('5/16/2009', '%m/%d/%Y')},
                       'Gonzo': {'subject_id': 'Gonzo', 'age': '9', 'description': 'NA',
                                 'genotype': 'NA', 'sex': 'M', 'species': 'Macaca mulatta', 'weight': '11',
                                 'date_of_birth': datetime.strptime('9/29/2008', '%m/%d/%Y')}}

    # File to import
    session_name = sessions.iloc[sess_iloc]['Session']
    subject = sessions.iloc[sess_iloc].Subject
    date = sessions.iloc[sess_iloc].Date
    fname = 'gatingwm_' + session_name.lower() + '.ns6'
    blk_file_name = path.join(fsroot, md.preproc_dat_superstruct, subject, date, fname)
    stim_path = path.join(fsroot, 'buschman', 'Projects', 'GatingInWorkingMemory', 'Files', 'Images', subject)


    #######################
    ## Initialize NWB file

    # Filter out current session
    units_df = units.loc[session_name]
    trials_df = trials.loc[session_name]
    conditions_df = conditions[conditions['Session'].eq(session_name)]
    events_df = events.loc[session_name]


    # Information about the experiment
    experiment = {'name': 'Gating In Working Memory',
                  'animal': subject,
                  'experimenter': 'Nikola Markov and Pavlos Kollias',
                  'lab': 'Buschman (www.timbuschman.com)',
                  'institution': 'Princeton Neuroscience Institute'}
    subject_info = subject_to_info[subject]
    animal_screen_distance = 61
    pixel_resolution_dva_per_pixel_x = 0
    pixel_resolution_dva_per_pixel_y = 0
    notes = 'Animal was positioned %d cm from the screen. Therefore, screen had resolution of %d and %d pixels per degree ' \
            'of visual angle along the horizontal and vertical axes, respectively.' % (animal_screen_distance,
                                                                                       pixel_resolution_dva_per_pixel_x,
                                                                                       pixel_resolution_dva_per_pixel_y)
    # Electrode groups should bring together of sites on the same electrode (i.e. for laminar arrays)
    electrode_group = [{'name': str(units_df.groupby('ChanGrp').get_group(grp)['ChanGrp'].values.unique()[0]),
                  'description': units_df.groupby('ChanGrp').get_group(grp)['ChanType'].values.unique()[0],
                  'location': units_df.groupby('ChanGrp').get_group(grp)['Area'].values.unique()[0]}
                 for grp
                 in units_df['ChanGrp'].unique()]

    # Information about electrodes -- can probably be extracted from ElectrodeTracker data
    electrode_info = [{'id': chan,
                       'impedance': 0.0,
                       'location': units_df.reset_index(drop=True).groupby('ChanNum').get_group(chan)['Area'].unique()[0],
                       'group': str(units_df.reset_index(drop=True).groupby('ChanNum').get_group(chan)['ChanGrp'].unique()[0])}
                      for chan
                      in units_df['ChanNum'].unique()]

    elec_list = units_df['ChanNum'].unique()

    # Initialize NWB file with information
    nwb_file_name = ConvertBlackrockToNWB.InitializeNWBFromBlackrock(blk_file_name, '',
                                                                     experiment=experiment,
                                                                     electrode_info=electrode_info,
                                                                     electrode_group=electrode_group,
                                                                     subject=subject_info,
                                                                     notes=notes,
                                                                     overwrite=True,
                                                                     verbose=True)


    ## Add raw stimulus images to file
    stimulus_image_info = [{'file': path.join(stim_path, '101.jpg'), 'name': 'Stimulus 1 (101)', 'resolution': None, 'description': 'First stimulus in Group 1'},
                           {'file': path.join(stim_path, '102.jpg'), 'name': 'Stimulus 2 (102)', 'resolution': None, 'description': 'Second stimulus in Group 1'},
                           {'file': path.join(stim_path, '201.jpg'), 'name': 'Stimulus 3 (201)', 'resolution': None, 'description': 'First stimulus in Group 2'},
                           {'file': path.join(stim_path, '202.jpg'), 'name': 'Stimulus 4 (202)', 'resolution': None, 'description': 'Second stimulus in Group 2'},
                           {'file': path.join(stim_path, '110.jpg'), 'name': 'Cue 1 (110)', 'resolution': None, 'description': 'Cue 1 for Group 1'},
                           {'file': path.join(stim_path, '120.jpg'), 'name': 'Cue 2 (120)', 'resolution': None, 'description': 'Cue 2 for Group 1'},
                           {'file': path.join(stim_path, '210.jpg'), 'name': 'Cue 3 (210)', 'resolution': None, 'description': 'Cue 1 for Group 2'},
                           {'file': path.join(stim_path, '220.jpg'), 'name': 'Cue 4 (220)', 'resolution': None, 'description': 'Cue 2 for Group 2'},
                           {'file': path.join(stim_path, '000.jpg'), 'name': 'Null Stimulus (000)', 'resolution': None, 'description': 'Behaviorally irrelevant stimulus'}]
    NWBHelperFunctions.AddStimulusImageToNWB(nwb_file_name, stimulus_image_info, verbose=True)



    #######################
    ## Add raw electrophysiology data
    ConvertBlackrockToNWB.AddBlackrockRawDataToNWB(blk_file_name, nwb_file_name, verbose=True, elec_ids=elec_list)


    #######################
    ## Add processed LFP data
    blk_lfp_file_name = path.splitext(blk_file_name)[0] + '.ns3'
    ConvertBlackrockToNWB.AddBlackrockLFPDataToNWB(blk_lfp_file_name, nwb_file_name, verbose=True, elec_ids=elec_list)


    #######################
    ## Filter LFP data from raw signal (not needed for most blackrock recordings)
    # NWBAnalysisFunctions.NWBFilterLFP(nwb_file_name, signal_name='Raw Data for Channels [1]', verbose=True, elec_ids=[1])


    #######################
    ## Add spiking data (from Plexon)

    # Get list of electrodes
    elec_list = units_df['ChanNum'].unique()

    # Create unit information dictionary
    unit_info = [{'electrode': unit_row['ChanNum'], 'unit': unit_row['UnitNum'],
                  'UserRating': unit_row['RatingCode'], 'UserComment': unit_row['RatingComment']}
                 for unit_ind, unit_row
                 in units_df.iterrows()]

    # Open file before adding spiking data to avoid bug in NWB format
    nwb_file, nwb_io = NWBHelperFunctions.OpenNWBFile(nwb_file_name)

    # Loop through electrodes, adding units from each Plexon file
    for cur_elec in elec_list:
        plx_file_name = path.join(fsroot, md.preproc_dat_superstruct, subject, date, 'CellSorting', '%s_chan%03d.plx' % (path.splitext(path.split(blk_file_name)[1])[0], cur_elec))
        ConvertPlexonToNWB.AddPlexonSpikeDataToNWB(plx_file_name, nwb_file=nwb_file, elec_ids=[cur_elec], add_units=True, unit_info=unit_info, verbose=True)

    # Write file and close
    NWBHelperFunctions.CloseNWBFile(nwb_file, nwb_io)


    #######################
    ## Add eye signals
    blk_cont_file_name = path.splitext(blk_file_name)[0] + '.ns4'
    ConvertBlackrockToNWB.AddBlackrockEyeDataToNWB(blk_cont_file_name, nwb_file_name, x_elec_id='ainp1', y_elec_id='ainp2', pupil_elec_id='ainp3', verbose=True) # TODO: Make sure the analog inputs are correct


    #######################
    # Add other analog signals (photodiode? bar position?)
    blk_cont_file_name = path.splitext(blk_file_name)[0] + '.ns4'
    signal_info = [{'label': 'ainp4', 'name': 'Photodiode', 'comments': 'Voltage output of photodiode monitoring screen changes.'}]
    ConvertBlackrockToNWB.AddBlackrockAnalogDataToNWB(blk_cont_file_name, nwb_file_name, signal_info=signal_info, verbose=True)  # TODO: Make sure the analog inputs are correct



    #######################
    ## Add trial information

    # Create key list for events and conditions
    # TODO: Check descriptions for the dataframe columns...
    # TODO: add indices pretty much to key_descriptions
    trial_key_description = {'TrialNum': 'The trial number in this session.',
                             'NumRewards': 'How many reward drops were given.',
                             'ReactionTime': 'Animal\'s reaction time, measured in Matlab as the time (in ms) from cue onset to response.',
                             'Cue': 'Cue identifier (= sub2ind([numGroups numCuesPerGroup],group,cue))',
                             'Stim': 'Stim identifier (= sub2ind([numGroups numStimPerGroup],group,stim))',
                             'FixAcq': 'TimeStamp of fixation acquisition.',
                             'FixBreak': 'TimeStamp of fixation break (if occurred).',
                             'BarAcq': 'Timestamp of bar holding.',
                             'BarRelease': 'Timestamp of bar release (if occurred).',
                             'RewardOn': 'TimeStamp (estimated) of reward period start.',
                             'RewardOff': 'TimeStamp (estimated) of reward period end.',
                             'StimOff': 'TimeStamp of last event display offset in trial',
                             'StageNum': 'How many stages (events) were in this trial.'}
    event_key_description = {'StageBitcodeIndex': 'TimeStamp of stage onset based on photodiode parsing',
                             'StageCategory': 'FixOnset, CueOnset, CueDelay, StimOnset, StimDelay, essentially unpacking of StageCode delay code, Nan for errors',
                             'StageCode': 'encode value coming from Matlab specs [2 11 10 8], NaN for errors',
                             'StageDuration': 'time until next event. Simply the difference between consecutive same trial event onset times',
                             'StageEncodeIndex': 'TimeStamp of stage onset based on encode triggers',
                             'StageIndex': 'Event number in trial event sequence, starts from 0 for fixation period',
                             'StageStim': 'Which cue/stimulus is on screen at current stage, -1 for fixation/delays',
                             'StageStimCategory': 'stimulus category - FIX for fixation and delay, CXX for cues, SXX for stimuli, NaN for errors',
                             'StageStimDelayCategory': 'CXX or SXX for stimulus preceding delay period, NaN for all other rows',
                             'StageTimeIndex': 'Official TimeStamp for stage onset to be used for analyses (either bitcode or extrapolated from encode). In ticks (so to convert to seconds, divide by sampling frequency)'}
    cond_key_description = {'StageIndex': 'Event number in trial event sequence, starts from 0 for fixation period',
                            'GatingCond': 'Event conditions label for Gating condition type. Levels: Cue, Gating, PreDist, PostDist, Target. Same label for stim and delay periods',
                            'GatingCondSpecialized': 'Specialized label from GatingCond: specializes CueDelay, GatingDelay, PreDistDelay, PostDistDelay, TargetDelay for delay periods',
                            'PostDistCategory': 'Event conditions label for type of postdistractor. Same if image is same group postdist, diff is image is diff group postdist, None for all other',
                            'StageStimExtended': 'Event conditions label for image (CXX or SXX) for current event',
                            'StageStimSpecialized': 'Specialized label from StageStimExtended: specializes CXXD and SXXD for delay periods'}

    # Iterate through trials and extract information
    trial_data = []
    for cur_trial_ind, cur_trial in trials_df.iterrows():
        # Take trial information from trials dataframe
        cur_trial_info = cur_trial[trial_key_description.keys()].to_dict()
        # Copy over the start and stop times
        cur_trial_info['start_time'] = cur_trial['StartTimeStamp']/30000
        cur_trial_info['stop_time'] = cur_trial['EndTimeStamp']/30000
        # Also create list of events during this trial
        trial_events = events_df.query('TrialNum==%d' % (cur_trial['TrialNum']))
        for cur_key in event_key_description.keys():
            cur_key_data = trial_events[cur_key].to_list()
            cur_trial_info[cur_key] = cur_key_data
        # Finally, create list of conditions during this trial
        trial_conds = conditions_df.query('TrialNum==%d' % (cur_trial['TrialNum']))
        for cur_key in cond_key_description.keys():
            cur_key_data = trial_conds[cur_key].to_list()
            cur_trial_info[cur_key] = cur_key_data
        # Append to overall list
        trial_data.append(cur_trial_info)

    # Add to NWB file
    NWBHelperFunctions.AddTrialDataToNWB(nwb_file_name, trial_data=trial_data, trial_descriptions={**trial_key_description, **event_key_description, **cond_key_description}, verbose=True)


    #######################
    ## Add stimulus epochs to file -- will track when something is on screen
    # Open NWB file for last additions
    nwb_file, nwb_io = NWBHelperFunctions.OpenNWBFile(nwb_file_name)

    # Define the sampling rate for the artificial signals
    sig_samp_freq = 30000

    # Combine the events and conditions dataframes into one UBERFRAME
    uber_df = pandas.merge(events_df.reset_index(drop=True), conditions_df.reset_index(drop=True))

    # The last stimulus doesn't have an offset, need to take it from the trials dataframe
    for cur_trial_ind, cur_trial in trials_df.iterrows():
        # Find the last stimulus, its duration is end of trial from trials_df minus the start
        trial_uber_ind = uber_df.index[uber_df['TrialNum'] == cur_trial['TrialNum']]
        if len(trial_uber_ind) == 0:
            continue
        if np.isnan(uber_df.at[trial_uber_ind[-1], 'StageDuration']):
            # Estimate from end of the stimulus sequence
            uber_df.at[trial_uber_ind[-1], 'StageDuration'] = cur_trial['StimOff'] - uber_df.at[trial_uber_ind[-1], 'StageTimeIndex']
        else:
            print("Oops, something wrong with parsing the uber dataframe.  Expecting a NaN at the end of each trial, but this didn't occur for trial %d." % (cur_trial['TrialNum']))


    # Define information about all of the stimulus signals to create
    # TODO: Check the names of these are correct -- specifically the Cue/Stim labels and their associations
    stim_info = [{'name': 'Fixation (before cue)', 'GatingCondSpecialized': [], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Cue 1 (Gate S11, S12)', 'GatingCondSpecialized': ['Cue'], 'StageStimCategory': ['C11'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Cue 2 (Gate S21, S22)', 'GatingCondSpecialized': ['Cue'], 'StageStimCategory': ['C21'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Cue 3 (Gate S11, S12)', 'GatingCondSpecialized': ['Cue'], 'StageStimCategory': ['C12'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Cue 4 (Gate S21, S22)', 'GatingCondSpecialized': ['Cue'], 'StageStimCategory': ['C22'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Cue 1 Delay', 'GatingCondSpecialized': ['CueDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory':['C11']},
                 {'name': 'Cue 2 Delay', 'GatingCondSpecialized': ['CueDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory':['C21']},
                 {'name': 'Cue 3 Delay', 'GatingCondSpecialized': ['CueDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory':['C12']},
                 {'name': 'Cue 4 Delay', 'GatingCondSpecialized': ['CueDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory':['C22']},
                 {'name': 'Stim 0', 'GatingCondSpecialized': ['PreDist', 'Gating', 'PostDist'], 'StageStimCategory': ['S00'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Stim 1', 'GatingCondSpecialized': ['PreDist', 'Gating', 'PostDist'], 'StageStimCategory': ['S11'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Stim 2', 'GatingCondSpecialized': ['PreDist', 'Gating', 'PostDist'], 'StageStimCategory': ['S21'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Stim 3', 'GatingCondSpecialized': ['PreDist', 'Gating', 'PostDist'], 'StageStimCategory': ['S12'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Stim 4', 'GatingCondSpecialized': ['PreDist', 'Gating', 'PostDist'], 'StageStimCategory': ['S22'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Stim 0 Delay', 'GatingCondSpecialized': ['PreDistDelay', 'GatingDelay', 'PostDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S00']},
                 {'name': 'Stim 1 Delay', 'GatingCondSpecialized': ['PreDistDelay', 'GatingDelay', 'PostDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S11']},
                 {'name': 'Stim 2 Delay', 'GatingCondSpecialized': ['PreDistDelay', 'GatingDelay', 'PostDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S21']},
                 {'name': 'Stim 3 Delay', 'GatingCondSpecialized': ['PreDistDelay', 'GatingDelay', 'PostDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S12']},
                 {'name': 'Stim 4 Delay', 'GatingCondSpecialized': ['PreDistDelay', 'GatingDelay', 'PostDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S22']},
                 {'name': 'Stim 0 PreDist', 'GatingCondSpecialized': ['PreDist'], 'StageStimCategory': ['S00'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Stim 1 PreDist', 'GatingCondSpecialized': ['PreDist'], 'StageStimCategory': ['S11'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Stim 2 PreDist', 'GatingCondSpecialized': ['PreDist'], 'StageStimCategory': ['S21'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Stim 3 PreDist', 'GatingCondSpecialized': ['PreDist'], 'StageStimCategory': ['S12'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Stim 4 PreDist', 'GatingCondSpecialized': ['PreDist'], 'StageStimCategory': ['S22'], 'StageStimDelayCategory':['NODELAY']},
                 {'name': 'Stim 0 PreDist Delay', 'GatingCondSpecialized': ['PreDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S00']},
                 {'name': 'Stim 1 PreDist Delay', 'GatingCondSpecialized': ['PreDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S11']},
                 {'name': 'Stim 2 PreDist Delay', 'GatingCondSpecialized': ['PreDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S21']},
                 {'name': 'Stim 3 PreDist Delay', 'GatingCondSpecialized': ['PreDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S12']},
                 {'name': 'Stim 4 PreDist Delay', 'GatingCondSpecialized': ['PreDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S22']},
                 {'name': 'Stim 0 Gated', 'GatingCondSpecialized': ['Gating'], 'StageStimCategory': ['S00'], 'StageStimDelayCategory': ['NODELAY']},
                 {'name': 'Stim 1 Gated', 'GatingCondSpecialized': ['Gating'], 'StageStimCategory': ['S11'], 'StageStimDelayCategory': ['NODELAY']},
                 {'name': 'Stim 2 Gated', 'GatingCondSpecialized': ['Gating'], 'StageStimCategory': ['S21'], 'StageStimDelayCategory': ['NODELAY']},
                 {'name': 'Stim 3 Gated', 'GatingCondSpecialized': ['Gating'], 'StageStimCategory': ['S12'], 'StageStimDelayCategory': ['NODELAY']},
                 {'name': 'Stim 4 Gated', 'GatingCondSpecialized': ['Gating'], 'StageStimCategory': ['S22'], 'StageStimDelayCategory': ['NODELAY']},
                 {'name': 'Stim 0 Gated Delay', 'GatingCondSpecialized': ['GatingDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S00']},
                 {'name': 'Stim 1 Gated Delay', 'GatingCondSpecialized': ['GatingDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S11']},
                 {'name': 'Stim 2 Gated Delay', 'GatingCondSpecialized': ['GatingDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S21']},
                 {'name': 'Stim 3 Gated Delay', 'GatingCondSpecialized': ['GatingDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S12']},
                 {'name': 'Stim 4 Gated Delay', 'GatingCondSpecialized': ['GatingDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S22']},
                 {'name': 'Stim 0 PostDist', 'GatingCondSpecialized': ['PostDist'], 'StageStimCategory': ['S00'], 'StageStimDelayCategory': ['NODELAY']},
                 {'name': 'Stim 1 PostDist', 'GatingCondSpecialized': ['PostDist'], 'StageStimCategory': ['S11'], 'StageStimDelayCategory': ['NODELAY']},
                 {'name': 'Stim 2 PostDist', 'GatingCondSpecialized': ['PostDist'], 'StageStimCategory': ['S21'], 'StageStimDelayCategory': ['NODELAY']},
                 {'name': 'Stim 3 PostDist', 'GatingCondSpecialized': ['PostDist'], 'StageStimCategory': ['S12'], 'StageStimDelayCategory': ['NODELAY']},
                 {'name': 'Stim 4 PostDist', 'GatingCondSpecialized': ['PostDist'], 'StageStimCategory': ['S22'], 'StageStimDelayCategory': ['NODELAY']},
                 {'name': 'Stim 0 PostDist Delay', 'GatingCondSpecialized': ['PostDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S00']},
                 {'name': 'Stim 1 PostDist Delay', 'GatingCondSpecialized': ['PostDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S11']},
                 {'name': 'Stim 2 PostDist Delay', 'GatingCondSpecialized': ['PostDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S21']},
                 {'name': 'Stim 3 PostDist Delay', 'GatingCondSpecialized': ['PostDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S12']},
                 {'name': 'Stim 4 PostDist Delay', 'GatingCondSpecialized': ['PostDistDelay'], 'StageStimCategory': ['FIX'], 'StageStimDelayCategory': ['S22']}]

    # Make sure behavior module exists in the file -- it really should
    if 'behavior' not in nwb_file.processing.keys():
        print("Adding behavior module...")
        signal_module = ProcessingModule(name='behavior', description="Processing module for behavior signals.")
        nwb_file.add_processing_module(signal_module)

    # Create a new behavioral events data interface
    stim_epoch_interface = BehavioralEpochs(name="Stimulus Image Epochs")

    # Loop through stimuli, creating sequences and adding them
    for cur_stim_ind, cur_stim_info in enumerate(stim_info):
        # Let's find the on time for each stimulus
        print("Finding onset and offset times for '%s'..." % (cur_stim_info['name']))
        query_str = ""
        key_list = np.setdiff1d([*cur_stim_info], ['name']).tolist()
        for cur_key in key_list:
            if len(cur_stim_info[cur_key]) > 0:
                if query_str != '':
                    query_str += ' and '
                query_str = query_str + "%s==%s" % (cur_key, cur_stim_info[cur_key])
        stim_events_df = uber_df.query(query_str)

        if stim_events_df.empty:
            warnings.warn("\tDidn't find any instances of '%s'. This may be expected, but throwing a warning anyways and skipping." % (cur_stim_info['name']), UserWarning)
            continue

        # Get onset and durations
        stim_events_time = stim_events_df[['StageTimeIndex', 'StageDuration']].to_numpy()
        # Convert durations to offsets
        stim_events_time[:, 1] = stim_events_time[:, 0] + stim_events_time[:, 1]

        # Create the onset/offset data matrix
        stim_events_data = np.zeros(stim_events_time.shape, dtype=int)
        stim_events_data[:, 0] = 1
        stim_events_data[:, 1] = -1

        # Flatten both time and data
        stim_events_time = stim_events_time.flatten()
        stim_events_data = stim_events_data.flatten()

        # Create interval series
        print("\tCreating interval time series...")
        cur_stim_series = IntervalSeries(name=cur_stim_info['name'],
                                         description='Epochs of when stimulus %s was on screen. This is extracted from trial information.' % (
                                         cur_stim_info['name']),
                                         data=stim_events_data,
                                         timestamps=stim_events_time)

        # Add to interface
        print("\tAdding signal to interface...")
        stim_epoch_interface.add_interval_series(cur_stim_series)

    # Add interface to the NWB file
    print("Adding stimulus epochs to behavior module.")
    nwb_file.processing['behavior'].add(stim_epoch_interface)

    # Write to file and close
    print("Writing NWB file and closing.")
    NWBHelperFunctions.CloseNWBFile(nwb_file, nwb_io)

    ###
    #######################
    ## Do wavelet transformation of the raw Add stimulus epochs to file -- will track when something is on screen
    # This will create a shallow copy of the file so the original remains untouched with analysis
    # NWBAnalysisFunctions.NWBWaveletTransform(nwb_file_name,
    #                                          signal_name=None,
    #                                          freq=np.concatenate((np.arange(1,5,1), np.arange(6,18,2), np.arange(18,40,4), np.arange(40,100,8))),
    #                                          width=2.5,
    #                                          elec_ids=None,
    #                                          wavelet_type='morlet',
    #                                          downsample=2,
    #                                          copy_file=True,
    #                                          verbose=True)


main()
