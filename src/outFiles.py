# Author: Katherina Cortes
# Date: June 14, 2022
# Purpose: To write information about viruses found in sample to csv file

# @params virusOutFile: file name and path to write to
# @params virusesFound: dict of virus IDs and counts
# @params viralHostsDB: dict of virus IDs and corresponding hosts
# @params viralNames: dict of virus IDs and name
# @writes: file containing virus ID, name, count and host information
#
# writes information about viruses found in sample to csv file
# formatted: virusID,virusName,count,host
def writeViruses(virusOutFile, virusesFound, viralHostsDB, viralNames):

    with open(virusOutFile, 'w') as of:
        for key in virusesFound.keys():
            if key in viralHostsDB:
                host = str(viralHostsDB[key])
                if ',' in host:
                    host = host.replace(',', ';')
            else:
                host = 'unknown'

            name = viralNames[key]
            if ',' in name:
                name = name.replace(',', ';')
            numViruses = virusesFound[key]
            of.write("%s,%s,%s,%s\n" % (key, name, numViruses, host))
    return