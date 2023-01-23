"""
20221206

make a config

input
    annotation (str) - split nullomers into +/- annotations using GENCODE annotations
    data_path (str) - full path to store data
    mutations (str) - full path to nullomers
    
method
    1. touch config file
    2. write key inputs
        - GENCODE annotation
        - DATA_PATH to store all data
        - NULLOMER data file. 
returns
"""

import os, sys
import configparser
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_/")

ANNOTATION = sys.argv[1] 
DATA_PATH = sys.argv[2]
MUTATIONS = sys.argv[3]
BUILD = sys.argv[4]
FLANK = sys.argv[5]

# import config reader/writer
import config_readwrite as crw


def touchConfig(annotation):
    
    """
    return config name, file to write to 
    
    inputs
        annotation (str) - for gencode
        
    method
        1. get path
        2. get config name, file name
        3. touch config file
    
    return 
        configname (str) - "config-<<type>>.ini"
        configfile (str) - FULL PATH + "config-<<type>>.ini"
        
    """
    #1
    path = os.getcwd()
    
    #2
    configname = "config-" + annotation +".ini"
    configfile = os.path.join(path, configname)
    
    #3
    cmd = "touch "+ configfile
    if os.path.exists(configfile) is False:
        os.system(cmd)
    
    return configname, configfile


def main(argv):
    
    # create the config file
    configname, configfile = touchConfig(ANNOTATION)
    
    # read the config file
    config, configfile_name = crw.read_config(configname)
    
    # add values to config file. 
    sections = ["DATA", "GENCODE", "RESULTS"]
    for s in sections:
        crw.check_section(config, s) 
    
    config["DATA"]["BUILD"] = BUILD
    config["DATA"]["PATH"] = DATA_PATH
    config["DATA"]["MUTATIONS"] = MUTATIONS
    config["DATA"]["FLANK"] = FLANK
    
        
    config["GENCODE"]["ANNOT"] = ANNOTATION
    config["GENCODE"]["VERSION"] = str(42)
    config["GENCODE"]["PATH"] = os.path.join(DATA_PATH, "Gencode")
    config["RESULTS"]["PATH"] = "/wynton/home/ahituv/fongsl/nullomers/results"
    
    
    # write the config file
    crw.write_config(config, configfile_name)

    print(configfile_name)
    
if __name__ == "__main__":
    main(sys.argv[1:])
