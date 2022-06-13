import os


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