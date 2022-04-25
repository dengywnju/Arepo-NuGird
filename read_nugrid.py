import h5py
import numpy as np



def readTable(filename):
    Table = {}
    Yields_data = []
    Ejected_Mass_Table = []
    LifeTime_Table = []
    Tflag = 0
    Masses = []
    Metallicities = []
    Species_names = []
    Yield_names = []
    Z0 = -1
    T = 0 #table flag
    yields_file = open(filename,'r')
    yields_lines = yields_file.readlines()
    for lines in yields_lines:
        if lines[0:7] ==  'H Table': 
            TableName = lines.split('=')
            M = float(TableName[1].split(',')[0])
            Z = float(TableName[2][:-2])
            Masses.append(M)
            if Z != Z0:
                Metallicities.append(Z)
                Yield_names.append('Z_'+str(Z))
                try:
                    Ejected_Mass_Table.append(np.array(Ejected_Mass_List))
                    LifeTime_Table.append(np.array(LifeTime_List))
                    Yields_data.append(Yields_Z)
                    LifeTime_List = []
                    Ejected_Mass_List = []
                    Yields_Z = []
                    Z0 = Z
                except:
                    LifeTime_List = []
                    Ejected_Mass_List = []
                    Yields_Z = []
                    Z0 = Z
                
        if lines[0:10] == 'H Lifetime': 
            LifeTime_List.append(float(lines.split(':')[1]))
        if lines[0:8] == 'H Mfinal':
            Mfinal = float(lines.split(':')[1])
            Ejected = M-Mfinal
            Ejected_Mass_List.append(Ejected)
        if lines[0] == '&':
            if lines[0:4] == '&Iso': 
                try:
                    Yields_Z.append(Yields_M_Z)
                    T = 1
                    Yields_M_Z = []
                except:
                    Yields_M_Z = []
                
            else:
                if T == 0:
                    Species_names.append(str(lines.split('&')[1]).strip())
                Yields = float(lines.split('&')[2])
                Yields_M_Z.append(Yields)
    Yields_Z.append(Yields_M_Z)
    Ejected_Mass_Table.append(np.array(Ejected_Mass_List))
    LifeTime_Table.append(np.array(LifeTime_List))
    Yields_data.append(Yields_Z)
    
    
    
    Table['Masses'] = np.array(Masses)
    Table['Metallicities'] = np.array(Metallicities)
    Table['Species_names'] = np.bytes_(Species_names)   
    Table['Yield_names'] = np.bytes_(Yield_names)
    Table['Yields'] = {}
    for i in range(len(Metallicities)):
        Table['Yields'][Yield_names[i]] = {}
        Table['Yields'][Yield_names[i]]['Yield'] = np.array(Yields_data[i]).T
        Table['Yields'][Yield_names[i]]['Ejected_Mass'] = np.array(Ejected_Mass_Table[i])
        Table['Yields'][Yield_names[i]]['LifeTime'] = np.array(LifeTime_Table[i])
    
    
    return Table


def findIsotope(speciesList,iso):
    return np.where(speciesList==np.bytes_(iso))[0][0]

### Rewrite original table!
def PopIIIturnToArepo(Table):    
    #turn to Arepo style y(m,Z)
    PrimordialAbundance = np.zeros(Table['Species_names'].shape[0])
    
    #Planck 18
    PrimordialAbundance[1] = 2*2.527e-5
    PrimordialAbundance[3] = 0.24672
    #Molaro 07
    PrimordialAbundance[2] = 1.05e-5*3
    PrimordialAbundance[5] = 7*4.41e-10

    PrimordialAbundance[0] = 1-np.sum(PrimordialAbundance[1:])
    returnTerm = np.multiply(PrimordialAbundance.reshape(PrimordialAbundance.shape[0],1),Table['Yields']['Z_0.0']['Ejected_Mass'].reshape((1,Table['Yields']['Z_0.0']['Ejected_Mass'].shape[0])),)
    
    Table['Yields']['Z_0.0']['Yield'] = Table['Yields']['Z_0.0']['Yield'] - returnTerm
    
    ind = findIsotope(Table['Species_names'],'Li-6') #find the beginning of metal
    Table['Yields']['Z_0.0']['Total_Metals'] = Table['Yields']['Z_0.0']['Yield'][ind:,:].sum(axis = 0)
    print(rawTable['Yields']['Z_0.0']['Total_Metals'])
    return Table

def PopIIIturnToArepoWith11Elements(Table): 
    
    #turn to Arepo style y(m,Z)
    PrimordialAbundance = np.zeros(Table['Species_names'].shape[0])
    
    #Planck 18
    PrimordialAbundance[1] = 2*2.527e-5
    PrimordialAbundance[3] = 0.24672
    #Molaro 07
    PrimordialAbundance[2] = 1.05e-5*3
    PrimordialAbundance[5] = 7*4.41e-10

    PrimordialAbundance[0] = 1-np.sum(PrimordialAbundance[1:])
    returnTerm = np.multiply(PrimordialAbundance.reshape(PrimordialAbundance.shape[0],1),Table['Yields']['Z_0.0']['Ejected_Mass'].reshape((1,Table['Yields']['Z_0.0']['Ejected_Mass'].shape[0])),)
    
    Table['Yields']['Z_0.0']['Yield'] = Table['Yields']['Z_0.0']['Yield'] - returnTerm
    
    ind = findIsotope(Table['Species_names'],'Li-6')#find the beginning of metal
    Table['Yields']['Z_0.0']['Total_Metals'] = Table['Yields']['Z_0.0']['Yield'][ind:,:].sum(axis = 0)
    
    # binning isotopes
    NewNames = np.bytes_(['Hydrogen', 'Helium', 'Carbon', 'Nitrogen', 'Oxygen', 'Neon',
       'Magnesium', 'Silicon', 'Sulphur', 'Calcium', 'Iron'])
    NewYileds = np.zeros((NewNames.shape[0],Table['Masses'].shape[0]))
    Hiso = ['H-1','H-2']
    Heiso = ['He-3','He-4']
    Ciso = ['C-12','C-13']
    Niso = ['N-14','N-15']
    Oiso = ['O-16','O-17', 'O-18']
    Neiso = ['Ne-20','Ne-21','Ne-22']
    Mgiso = ['Mg-24','Mg-25','Mg-26']
    Siiso = ['Si-28','Si-29','Si-30']
    Siso = ['S-32','S-33','S-34','S-36']
    Caiso = ['Ca-40', 'Ca-42', 'Ca-43', 'Ca-44', 'Ca-46', 'Ca-48']
    Feiso = ['Fe-54', 'Fe-56', 'Fe-57', 'Fe-58']
    Ele = [Hiso, Heiso, Ciso, Niso, Oiso, Neiso, Mgiso, Siiso, Siso, Caiso, Feiso]
    
    lines = 0
    for spec in Ele:
        for isotope in spec:
            where = findIsotope(Table['Species_names'],isotope)
            NewYileds[lines]+=Table['Yields']['Z_0.0']['Yield'][where,:]
        lines += 1
    Table['Species_names'] = NewNames
    Table['Yields']['Z_0.0']['Yield'] =  NewYileds
    return Table
        