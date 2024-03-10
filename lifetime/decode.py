#########################################################
#                                                       #
#                                                       #
#    Haider Ali                                         #
#    Tissue Mechanics and Remodeling Laboratory         #
#    Jacopo Ferruzzi                                    #
#    The University of Texas at Dallas                  #
#                                                       #
#                                                       #
#    A few functions to decode ISS fbd file             #  
#                                                       #
#                                                       #
#                                                       #
#                                                       #
#                                                       #
#     last edit: 08/21/2023                             #
#########################################################













import numpy as np
from numpy.typing import NDArray
import cv2
import xml.etree.ElementTree as ET
import os
import time
import matplotlib.pyplot as plt
from matplotlib.pyplot import imshow
import re
from tqdm import tqdm, trange
from dataclasses import dataclass


@dataclass
class Params:
    scan_line_left_border: int
    """the left border of the scan line."""
    scan_line_right_border: int
    """the right border of the scan line."""
    scan_line_length: int
    """the length of the scan line."""
    xpixels: int
    """the number of pixels in the x-axis."""
    ypixels: int
    """the number of pixels in the y-axis."""
    zpixels: int
    """the number of pixels in the z-axis."""
    dwell_time: float
    """the pixel dwell time in seconds."""
    excit_freq: int
    """the excitation frequency in Hz."""
    firmware: str
    """the path to the firmware file used in the session."""
    decoder: str
    """the decoder used in the session."""
    second_harmonic: bool
    """whether the second harmonic is used."""
    channel_mapping: str
    """the mapping of channels used in the session."""
    desc: str
    unit_time: float







@dataclass
class Result:
    phase: NDArray
    frame_signal: NDArray
    photons: NDArray[np.uint8]
    photon_win: NDArray[np.int8]
    micro: NDArray[np.uint8]
    macro: NDArray[np.uint64]
    macro_base: NDArray[np.uint64]
    elapsed_time_s: NDArray[np.float64]
    frame_time_s: NDArray[np.float64]






def extract_meta(acq_info_header_bytes):

    """
        returns
        - params : Params

        Takes in the header bytes and uses XML tree get the aquisition information
    """
    
    acq_info = acq_info_header_bytes.strip(b'\x00').decode('UTF-8')

    # XML element
    root = ET.fromstring(acq_info)

    # Find the ScanParams element
    scan_params_elem = root.find(".//ScanParams")
    firm_params_elem = root.find(".//FirmwareParams")

    # Extract the image parameters - add other params here if you wished
    params = Params(
        scan_line_left_border = int(scan_params_elem.find("ScanLineLeftBorder").text),
        scan_line_right_border = int(scan_params_elem.find("ScanLineRightBorder").text),
        scan_line_length = int(scan_params_elem.find("ScanLineLength").text),
        xpixels = int(scan_params_elem.find("XPixels").text),
        ypixels = int(scan_params_elem.find("YPixels").text),
        zpixels = int(scan_params_elem.find("ZPixels").text),
        dwell_time = float(scan_params_elem.find("PixelDwellTime").text),
        excit_freq = int(scan_params_elem.find("ExcitationFrequency").text),
        firmware = (firm_params_elem.find("FileName").text),
        decoder = (firm_params_elem.find("DecoderName").text),
        second_harmonic = bool(firm_params_elem.find("Use2ndHarmonic").text),
        channel_mapping = (firm_params_elem.find("ChannelMapping").text),
        desc = (firm_params_elem.find("Description").text),
        unit_time = -1,
    )

    params.desc = params.desc[:params.desc.find('/Vers')].split('/')

    if scan_params_elem.find("PixelDwellTime").attrib['Unit'] == 'Millisecond':
        params.dwell_time *= .001

    return params















def build_int_singlez(params, result, sec, channel=1):

    """
        Builds intensity images for a single frame (z slice). Mostly just used initially before we started building the decays.
    """
    
    frame_time = (result.frame_time_s[sec[0]:sec[1], channel] if sec else result.frame_time_s[:, channel])
    photons = (result.photons[sec[0]:sec[1], channel] if sec else result.photons[:, channel])

    # calculate how many pixels have been imaged at each time point in frame time 
    n_pix = (frame_time // params.dwell_time)

    # Calculate the corresponding row/column of each pixel number
    row = (n_pix // params.scan_line_length).astype(int)
    col = ((n_pix % params.scan_line_length) - params.scan_line_left_border).astype(int)
    # col = (n_pix % params.scan_line_length).astype(int)
    
    # find all the photons that arrived and their pixel
    valid_photons = (
        (col >= 0) & (col < params.xpixels) & (row < params.ypixels) & photons
    ).astype(bool)

    image = np.zeros((params.ypixels, params.xpixels), dtype=np.uint64)

    # use np.add.at to accumulate values into the image array. this is a vectorized solution to the for loop beneath it
    np.add.at(image, (row[valid_photons], col[valid_photons]), 1)
    # for i in trange(len(valid_photons), desc="build_int_single_z"):
    #     if valid_photons[i]:
    #         image[row[i], col[i]] += 1


    return image











def _build_images_seperate_avg(params, result, channel=1):
    
    """
        Same as the build images singlez but this will seperate the frames that we have averaged over. Primarily was just used for debugging. 
    """

    image = np.zeros((params.ypixels, params.xpixels, np.sum(result.frame_signal), params.zpixels), dtype=np.uint16) 
    frame_sigs = np.where(result.frame_signal == 1)[0]
    
    for z in range(params.zpixels):

        for frame in range(len(frame_sigs)-1):
            slice_start = frame_sigs[frame]
            slice_end = frame_sigs[frame+1] - 1
            frame_time = result.frame_time_s[slice_start:slice_end, channel]
            image[:,:,frame,z] = build_int_singlez(frame_time, result.photons[slice_start:slice_end, channel], params)
         
    return image
       











def build_decays(params, result, channel=1):

    """
        This function takes in all the decoded data and then seperates the total experiment time based on which frame is being built. 
        The frame time and micro time are used to then build decays for each pixel by calling the build_decays_singlez. 
    """

    
    # image = np.zeros((params.xpixels, params.ypixels, params.zpixels), dtype=np.uint16)
    decays = np.zeros((params.xpixels, params.ypixels, params.zpixels, 128), dtype=np.uint16)

    # subtract frame per z by 1 to get rid of last layer (avg)
    # !!!!!!! need to check if the last layer has reached the end of the image
    frame_per_z = int(np.sum(result.frame_signal) / params.zpixels)
    frame_sigs = np.where(result.frame_signal == 1)[0]
    frame_sigs = np.concatenate((frame_sigs, [result.frame_signal.shape[0]]))
    
    for z in trange(params.zpixels, desc="z"):
        
        slice_start = frame_sigs[z * frame_per_z]
        slice_end = frame_sigs[((z+1) * frame_per_z)] - 1
        # image[:,:,z] = build_int_singlez(params, result, sec=[slice_start, slice_end])
        decays[:,:,z,:] = build_decays_singlez(params, result, sec=[slice_start, slice_end])

    return decays




















def build_decays_singlez(params, result, channel=1, sec=[]):

    """
        This function uses the frame time and other params to build decays. It takes in 5 parameters:

            - frame_time: this is the frame times for all the frames averaged for a single z slice (16 frames for any LSM image or MATL)
            - result: 
            - params: 
            - channel: which channel we are decoding (default is 1 because mostly use NADH)
            - sec: these are the data start and end indicies of the slice of interest 
    """
    
    frame_time = (result.frame_time_s[sec[0]:sec[1], channel] if sec else result.frame_time_s[:, channel])
    photons = (result.photons[sec[0]:sec[1], channel] if sec else result.photons[:, channel])

    n_pix = (frame_time // params.dwell_time)

    row = (n_pix // params.scan_line_length).astype(int)
    col = ((n_pix % params.scan_line_length) - params.scan_line_left_border).astype(int)

    bphotons = photons == 1  
    valid_photons = (
        (col >= 0) & (col < params.xpixels) & (row < params.ypixels) & bphotons
    ).astype(bool)
    micro = result.micro[sec[0]:sec[1], channel] 

    decays = np.zeros((params.xpixels, params.ypixels, 128), dtype=np.uint64)

    # vectorized approach to the for loop beneath it
    np.add.at(decays, (row[valid_photons], col[valid_photons], micro[valid_photons]), 1)
    # for i in trange(len(valid_photons), desc="build_decays_singlez"):
    #     if valid_photons[i]:
    #         decays[row[i], col[i], micro[i]] += 1

    return decays

































def decode_8w(data, params):


    """
        This function is used to decode the fbd file and read the file header.

        Parameters:
        - path: the path to the fbd file

        Returns:
        - result (dict)
            go look at the definition

    """

    # Define classical decoder
    classical_decoder = np.array([
        0x00,
        0x01, 0x03, 0x05, 0x07, 0x09, 0x0b, 0x0d, 0x0f,  # w0-w7, ch1
        0x10, 0x30, 0x50, 0x70, 0x90, 0xb0, 0xd0, 0xf0,  # w0-w7, ch2
        0x11, 0x13, 0x15, 0x17, 0x19, 0x1b, 0x1d, 0x1f,  # w0-w7, ch1, w0 ch2
        0x31, 0x33, 0x35, 0x37, 0x39, 0x3b, 0x3d, 0x3f,  # w0-w7, ch1, w1 ch2
        0x51, 0x53, 0x55, 0x57, 0x59, 0x5b, 0x5d, 0x5f,  # w0-w7, ch1, w2 ch2
        0x71, 0x73, 0x75, 0x77, 0x79, 0x7b, 0x7d, 0x7f,  # w0-w7, ch1, w3 ch2
        0x91, 0x93, 0x95, 0x97, 0x99, 0x9b, 0x9d, 0x9f,  # w0-w7, ch1, w4 ch2
        0xb1, 0xb3, 0xb5, 0xb7, 0xb9, 0xbb, 0xbd, 0xbf,  # w0-w7, ch1, w5 ch2
        0xd1, 0xd3, 0xd5, 0xd7, 0xd9, 0xdb, 0xdd, 0xdf,  # w0-w7, ch1, w6 ch2
        0xf1, 0xf3, 0xf5, 0xf7, 0xf9, 0xfb, 0xfd, 0xff   # w0-w7, ch1, w7 ch2
    ], dtype=np.uint8)

    # extract phase, frame signal, and compress code, then uncompress
    phase = data & 0xFF
    frame_signal = np.bitwise_and(data >> 8, 0x01)
    compress_code = np.bitwise_and(data >> 9, 0x7F)
    
    decompress = classical_decoder[compress_code].astype(np.uint8)

    # extract photon information from the decoded data
    micro_max = 2 ** 8
    channel_count = 2
    params.unit_time = (micro_max - 1) / (micro_max * params.excit_freq) * (2 if params.second_harmonic  else 1) # time of 256 bins in s
    phase_zero = phase == 0
    macro_time_base = np.zeros_like(phase, dtype=np.uint64)
    macro_time_base = np.cumsum(np.where(phase_zero, micro_max, 0)) # the np.where evaluates the condition phase_zero and asigns 256 to true and 0 otherwise
    frame_indicies = np.where(frame_signal == 1)[0]
    phase_bins = micro_max 
    window_bins = 2 ** 3
    display_bins = phase_bins / (2 if params.second_harmonic else 1) 

    # initialize result
    result = Result(
        phase = phase,
        frame_signal = frame_signal,
        photons = np.zeros((data.size, channel_count), dtype=np.uint8),
        photon_win = np.full((data.size, channel_count), -1, dtype=np.int8),
        micro = np.full((data.size, channel_count), -1, dtype=np.uint8),
        macro = np.full((data.size, channel_count), -1, dtype=np.uint64),
        macro_base = np.zeros((data.size,), dtype=np.uint64),
        elapsed_time_s = np.zeros((data.size, channel_count), dtype=np.float64),
        frame_time_s = np.zeros((data.size, channel_count), dtype=np.float64),
    )


    for channel in range(channel_count):

        # 
        result.photons[:, channel] = np.bitwise_and((decompress >> (4 * channel)), 0x01)
        result.photon_win[:, channel] = np.bitwise_and((decompress >> (1 + 4 * channel)), 0x07)
        
        # result.micro[:, channel] = ((display_bins - 1 - (phase + (result.photon_win[:, channel] * phase_bins/window_bins))) % display_bins) // (2 if params.second_harmonic else 1)  
        result.micro[:, channel] = ((display_bins - 1 - (phase + (result.photon_win[:, channel] * phase_bins/window_bins))) % display_bins)
        
        # result.macro[:, channel] = macro_time_base + result.micro[:, channel] # going to keep this line because it was the only thing wrong for ages
        result.macro[:, channel] = macro_time_base + phase

        # experiment time in seconds
        result.elapsed_time_s[:, channel] = result.macro[:, channel] * params.unit_time

        # seperate total expriment time into frame times
        for i in range(len(frame_indicies) - 1):
            
            start = frame_indicies[i]
            end = frame_indicies[i+1] - 1

            frame_start = result.elapsed_time_s[start, channel]
            result.frame_time_s[start:end, channel] = (result.elapsed_time_s[start:end, channel] - frame_start)
        
        frame_start = result.elapsed_time_s[(end+1), channel]
        result.frame_time_s[(end+1):, channel] = (result.elapsed_time_s[(end+1):, channel] - frame_start)


    result.macro_base = macro_time_base


    return (params, result)














def main(path):

    """
        This is the function that matlab will call. It takes in the path to the fbd file and decodes it and returns aquisition information and the decays.

        Current supported decoders:
            - Digilent_2Ch_8W_Div2 (8w)
    """

    tic1 = time.perf_counter()
    if '$XX2X' in path:

        # Open the file and read the header
        with open(path, "rb") as f:

            # Read 1KB type header
            f.read(1024)
            # Read 32KB Acq Info Header with info starting at 16KBytes/32KB Location
            acq_info_header_bytes = f.read((32*1024))

            # Get the size of the data section
            data_size_bytes = os.path.getsize(path) - (33*1024)

            # Read the data section (> is  big endian)
            data = np.fromfile(f, dtype='<u2', count=data_size_bytes)
        
        params = extract_meta(acq_info_header_bytes)
    
    else:
        raise Exception('Unknown number of header data bytes')
    
    
    if params.decoder != '8w':
        raise Exception('Decoder is not 8w')
    else:
        params, result = decode_8w(data, params)
        decays = build_decays(params, result) 
        # image = build_images_seperate_avg(result, params)


    toc1 = time.perf_counter()
    print(f"took: {toc1 - tic1:0.4f} seconds")

    return params, decays
















if __name__ == "__main__":
    
    # path = '/Users/haider/Documents/tmr/data/20230203/ISS/calibration/Fluoresscein_Cali_740nm_750hv$XX2X.fbd'
    path = '/Users/haider/Documents/tmr/data/20230203/ISS/control/monolayer2/area1/740nm/MCF10A_CNTRL_-D_control_monolayer2_area1_740nm$XX2X.fbd'
    params, decays = main(path)

