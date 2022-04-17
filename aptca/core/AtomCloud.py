"""
Atom point cloud for APT data Analysis
"""
from enum import Enum
import os
import struct
import sys
from typing import Any, Tuple, Union
from warnings import warn
import warnings
import numpy as np
import pandas as pd
from sklearn.neighbors import KDTree
from scipy.stats import mode
from aptca.core.RRNG import RRNG

from aptca.utils import validate
from aptca.utils.logging import log
from aptca.utils import quickplots



class AtomCloud:
    """
    3D point cloud data structure based on K-d tree
    """
    def __init__(self):
        """
        Initialize point cloud class. Taking a numpy array m by n size, \
        where m is the number of data points, n is the attribute of each data point.

        Note data size is not checked, it can be anything eligible.
        """
        self._pos = None
        self._mass = None
        self._mass_random = None
        self._num_data = None
        self._data_dim = None
        self._loaded = False
        self._name = None
        self._path = None
        
    def load_atomcloud(self, pos, mass):
        if not isinstance(pos, np.ndarray) or not isinstance(mass, np.ndarray):
            raise TypeError("Mass and xyz coordinates must be numpy arrays")
        if len(pos.shape) != 2 or pos.shape[1] != 3:
            raise ValueError(f"XYZ array is not correct shape {pos.shape} should be (Mx3)")
        if len(mass.shape) != 1:
            print("Mass shape", len(mass.shape))
            raise ValueError("Mass array must be 1 dimensional")
        if pos.shape[0] != mass.shape[0]:
            raise ValueError("XYZ and mass must have the same number of entries (per atom)")

        self._pos = pos
        self._mass = mass
        self._mass_random = np.random.permutation(mass)
        self._num_data, self._data_dim = self._pos.shape
        self.gen_kdtree()
        self._loaded = True
    
    def load_atomcloud_fromfile(self,apt_file:str):
        validate.file_exists(apt_file)
        self._name = '.'.join(os.path.basename(apt_file).split('.')[0:-1])
        self._path = os.path.dirname(apt_file)
        extension_atom = apt_file.split('.')[-1]
        if extension_atom == 'POS' or extension_atom == 'pos':
            pos,mass = AtomCloud.from_pos(apt_file)
            self.load_atomcloud(pos,mass)
            self._loaded = True
        elif extension_atom == 'APT' or extension_atom == "apt":
            pos,mass = AtomCloud.pos_from_apt(apt_file)
            self.load_atomcloud(pos,mass)
            self._loaded = True
        else:
            log.error(f"File format not supported: {extension_atom}. Choose from .pos/.POS or .apt/.APT file")
    
    def map_mass_to_ion(self,rrng:RRNG):
        # slow (super slow)
        self._ion_from_mass = np.array(list(map(rrng.get_ion_from_mass,self._mass)))

    def gen_atom_density_distribution(self,plot=False,save=True,save_path=None):
        if not self._loaded:
            raise ValueError("Please load AtomClooud first either from file or from pos,mass data!")
        
        x = self.pos[:,0]
        y = self.pos[:,1]
        z = self.pos[:,2]
        # move to positive coord system
        X = np.round(x-min(x)).astype(int)
        Y = np.round(y-min(y)).astype(int)
        Z = np.round(z-min(z)).astype(int)
        # ranges for positive integer lists of X Y Z 
        X_range = np.linspace(min(X),max(X),num=max(X)-min(X)+1,dtype=int)
        Y_range = np.linspace(min(Y),max(Y),num=max(Y)-min(Y)+1,dtype=int)
        Z_range = np.linspace(min(Z),max(Z),num=max(Z)-min(Z)+1,dtype=int)
        n_atoms_at_location = np.zeros((len(X_range),len(Y_range),len(Z_range))).astype(int)
        # count atoms in locations [i-0.5,i+0.5), i is integer and from X, Y, or Z.
        for i in range(len(X)):
            n_atoms_at_location[X[i],Y[i],Z[i]] += 1
        
        #return n_atoms_at_location
        
        D = np.zeros(len(X_range)*len(Y_range)*len(Z_range)).astype(int)
        Idx = np.zeros(len(X_range)*len(Y_range)*len(Z_range)).astype(int)
        Idy = np.zeros(len(X_range)*len(Y_range)*len(Z_range)).astype(int)
        Idz = np.zeros(len(X_range)*len(Y_range)*len(Z_range)).astype(int)
        index = 0
        counter = 0
        for i in range(1,max(X_range)-1):
            for  j in range(1,max(Y_range)-1):
                for k in range(1,max(Z_range)-1):
                    # Determine if surrounding voxels are empty
                    for n in range(i-1,i+2):
                        for q in range(j-1,j+2):
                            for r in range(k-1,k+2):
                                if n_atoms_at_location[n,q,r] == 0:
                                    counter += 1
                    #If no surrounding voxels are not empty, add this voxel to the list
                    if counter == 0:
                        Idx[index] = i
                        Idy[index] = j
                        Idz[index] = k
                        D[index] = n_atoms_at_location[i,j,k]
                        index = index+1
                    # reset counter for surrounding atoms counting
                    counter =0
                    
        # Remove empty voxels from the list
        non_empty_voxels = np.where(D > 0)
        #print(non_empty_voxels)
        Density = D[non_empty_voxels]
        Idx = Idx[non_empty_voxels]
        Idy = Idy[non_empty_voxels]
        Idz = Idz[non_empty_voxels]
        # Calculation of Density vs distance_axial and Density vs Distance from central Axis
        # Creating Density plot with distance_axial
        #return D, Density, Idx, Idy, Idz, non_empty_voxels, n_atoms_at_location
        distance_axial = np.zeros(max(Idz)).astype(int)
        avg_density_axial = np.zeros(max(Idz))
        for i in range(max(Idz)):
            distance_axial[i] = i
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                avg_density_axial[i] = np.mean(Density[np.where(Idz == i)])/0.52
        # Creating Density plot with radial
        # Find Central Axis
        AxisX = (max(x) + min(x))/2-min(x)
        AxisY = (max(y) + min(y))/2-min(y)
        Rdistance = np.sqrt(np.power(Idx - AxisX,2)+np.power(Idy - AxisY,2))

        distance_radial_all = np.round(Rdistance).astype(int)
        distance_radial = np.linspace(0,max(distance_radial_all)-1,num=max(distance_radial_all),dtype=int)
        avg_density_radial = np.zeros(max(distance_radial_all))
        for i in distance_radial:
            avg_density_radial[i] = np.mean(Density[np.where(distance_radial_all == i)])/0.52
        
        calc_density = mode(Density)[0]/0.52
        print(f"Calc_density = {calc_density}")

        length_axial = len(avg_density_axial)
        axial10Percent = round(0.1*length_axial)
        axial80Percent = round(0.8*length_axial)
        axial_density_Range_1_8 = avg_density_axial[axial80Percent] - avg_density_axial[axial10Percent]
        axial_range = max(avg_density_axial[axial10Percent:axial80Percent]) - min(avg_density_axial[axial10Percent:axial80Percent])
        length_radial = len(avg_density_radial)
        radial10Percent = round(0.1*length_radial)
        radial80Percent = round(0.8*length_radial)
        radial_density_Range_1_8 = avg_density_radial[radial80Percent] - avg_density_radial[radial10Percent]
        radial_range = max(avg_density_radial[radial10Percent:radial80Percent]) - min(avg_density_radial[radial10Percent:radial80Percent])
        print(f"Axial density range: {axial_density_Range_1_8}")
        print(f"Radial density range: {radial_density_Range_1_8}")
        print("correct ranges:")
        print(f"Axial density range: {axial_range}")
        print(f"Radial density range: {radial_range}")
        self.result = n_atoms_at_location
        self.D = D
        self.Idx = Idx
        self.Idy = Idy
        self.Idz = Idz
        self.Rdistance = Rdistance
        self.Density = Density
        self.distance_axial = distance_axial
        self.avg_density_axial = avg_density_axial
        self.distance_radial = distance_radial
        self.avg_density_radial = avg_density_radial
        self.axial_range = axial_range
        self.radial_range = radial_range
        self.axial_density_Range_1_8 = axial_density_Range_1_8
        self.radial_density_Range_1_8 = radial_density_Range_1_8
        if save:
            if save_path == None:
                save_path = self._path
        save_name = self._name + '-DensityDistribution'
        if plot:
            quickplots.plot_density_distribution(Density,distance_axial,avg_density_axial,
                                                distance_radial,avg_density_radial,save=save,save_path=save_path,save_name=save_name)
        #return Density, distance_axial, avg_density_axial, distance_radial, avg_density_radial

        
    @property
    def pos(self) -> np.ndarray:
        return self._pos
    
    @property
    def mass(self) -> np.ndarray:
        return self._mass

    @property
    def mass_random(self) -> np.ndarray:
        return self._mass_random

    @property
    def total_atoms(self):
        """
        Return total number of data points .
        """
        return self._num_data

    def get_coord(self, idx):
        """
        Return coordinates of point with index idx.
        """
        return self._pos[idx]

    def gen_kdtree(self, metric='euclidean'):
        """
        Generate KDTree for datapoints. Using SciPy.neighbors.KDTree method.
        Default distance metric is set to euclidean.
        """
        self.kdtree = KDTree(self._pos, metric=metric)

    def query_neighbors(self, idx, minpts):
        """
        Obtain distance and index information for minpts-neighbors of the point with index idx.
        It is bassically the same as sklear.neighbors.kdtree.query, but instead \
        of return a 1-by-n numpy array, the results are reshape to simple 1D numpy \
        array of size minpts.

        idx:
            integer, index for the point of interest in point cloud.
        minpts:
            integer, the number of nearest neighbor to return. Note that for our case k=1 return distance and index \
            of itself, since all our query points are in the point cloud.
        Return:
            (dist, ind), each is an 1D numpy array of size minpts.
        """
        dist, ind = self.kdtree.query(self._pos[idx].reshape(-1, 3), k=minpts)
        # note the dist, ind could be multi-dimensional array, with each row correspond to a sample queried. \
        # in this function, only ONE sample is allowed for query, so that return dist[0] to reduce dimension to 1D.
        return dist[0], ind[0]

    def query_radius(self, idx, eps):
        """
        Return distance and index information for neighbors within eps distance to the point with idx .

        idx:
            integer, index for the point of interest in point cloud.
        eps:
            float, the search distance to find nearest neighbors.
        Return:
            (dist, ind), each is an 1D numpy array of size minpts.
        """
        # note the dist, ind could be multi-dimensional array, with each row correspond to a sample queried. \
        # in this function, only ONE sample is allowed for query, so that return dist[0] to reduce dimension to 1D.
        ind, dist = self.kdtree.query_radius(self._pos[idx].reshape(-1, 3), r=eps, return_distance=True)
        return dist[0], ind[0]

    def get_KNN_dist(self, num_bins, k):
        """
        Obtain KNN distribution. This is only a wrapper over query_neighbor method and numpy.histogram function.
        """
        idx = self.total_atoms-1
        dist, _ = self.query_neighbors(idx, k)
        hist, bin_edge = np.histogram(dist, num_bins) #Note bin_edge is len(hist) + 1
        return hist, bin_edge


    @classmethod
    def pos_from_txt(cls,filename, sep="\t", header=None):
        """
        Read txt data file. Assuming the first line is file header and the 
        rest m row x n col are data. There should be no missing values. Delimiter is
        assumed to be space.

        It basically calls numpy.loadtx set delimiter default to space and skiprows default to 1
        for my laziness to input those two parameters everytime.
        
        filename:
            filename of the file want to be read. Default data is float type.
        header:
            number of rows that considered as header of the file.
        """
        return pd.read_csv(filename, sep=sep,header=header).to_numpy()
   
    @classmethod
    def write_pos_to_txt(cls,filename, data, header='X Y Z Da'):
        """
        Write data to txt format. Basically calls np.savetxt, again just a minor helper to keep all output consistent format.
        
        data:
            a numpy array of m-by-n size
        filename:
            a string for filename to store data
        header:
            a string serves as file header   
        """
        np.savetxt(filename, data, delimiter=' ', header=header, comments='')
        return

    @classmethod
    def from_pos(cls,fname):
        """
        Loads an APT .pos file and return a column-based data numpy array.

        Output formate is:
        [[x0, y0, z0, Da0],
        [x1, y1, z1, Da1],
        ...]
        """
        dt = np.dtype('>f4') #default pos file format is big endian byteorder four byte float. (check this in reference book)
        d = np.fromfile(fname, dtype=dt, count=-1) #use numpy from file to read a binary file. The d format is [x0,y0,z0,D0,x1,y1,z1,D1,...]
        data = np.reshape(d, (-1, 4)) # change the data shape to [m rows, 4 cols] from a total data size of 4m.
        pos = data[:,0:3]
        mass = data[:,3]
        #aptdata =cls(pos,mass)
        return pos,mass

    @classmethod
    def write_pos(cls,pos_fname, data):
        """
        Writing a numpy array data to APT readable .pos format.

        data:
            a n by 4 numpy array [[x0,y0,z0,Da0],
                                [x1,y1,z1,Da1],
                                ...]
        pos_fname:
            filename for the output pos.
        """
        if data.shape[1] != 4:
            sys.exit('data must be a numpy array of size m-by-4')

        assert data.dtype == '>f4' # very important, default float datatype will give corrupted result.

        flat_data = np.ndarray.flatten(data)
        flat_data.tofile(pos_fname) # Note to self, format in tofile dose not work for changing data type. 
                                    # So it need an assertation ahead to make sure data type is corrects.
        
        print('Pos file writing finished')   
        return

    
    @classmethod
    def pos_from_apt(cls, filepath: str, verbose: bool = False):
        """
        Read the contents of an apt file into a Roi container

        :param filepath: Path to apt file
        :param verbose: Print the structure of the apt file as it is read (for debug purposes)
        """

        validate.file_exists(filepath)
        log.info("Reading apt file: {}".format(os.path.basename(filepath)))

        class RelType(Enum):
            REL_UNKNOWN = 0
            ONE_TO_ONE = 1
            INDEXED = 2,
            UNRELATED = 3
            ONE_TO_MANY = 4

        class RecordType(Enum):
            RT_UNKNOWN = 0
            FIXED_SIZE = 1
            VARIABLE_SIZE = 2
            VARIABLE_INDEXED = 3

        class RecordDataType(Enum):
            DT_UNKNOWN = 0
            INT = 1
            UINT = 2
            FLOAT = 3
            CHARSTRING = 4
            OTHER = 5

        class Dtype(Enum):
            int32 = 4
            int64 = 8
            char = 1
            wchar_t = 2
            filetime = 8

        def record_dtype2numpy_dtype(rec_dtype: RecordDataType, size: int):
            """
            Map a section's record data type to its equivalent numpy dtype
            """
            if rec_dtype in (RecordDataType.UINT, RecordDataType.CHARSTRING):
                raise ValueError("Cannot map to UINT or CHARSTRING")

            int_map = {
                8: np.int8,
                16: np.int16,
                32: np.int32,
                64: np.int64
            }

            float_map = {
                32: np.float32,
                64: np.float64
            }

            if rec_dtype == RecordDataType.INT:
                return int_map[size]
            elif rec_dtype == RecordDataType.FLOAT:
                return float_map[size]
            else:
                raise ValueError(f"Unexpected record data type {rec_dtype}")

        # Maps the apt format data type to str format needed for struct.unpack
        dtype2fmt = {
            Dtype.int32: "i",
            Dtype.int64: "q",
            Dtype.char: "c",
            Dtype.filetime: "Q",
            Dtype.wchar_t: "c"
        }

        # Maps the apt format data type to python data type
        dtype2constr = {
            Dtype.int32: int,
            Dtype.int64: int,
            Dtype.char: lambda x: x.decode("utf-8"),
            Dtype.wchar_t: lambda x: x.decode("utf-16"),
            Dtype.filetime: int
        }

        with open(filepath, "rb") as dat:
            def read_chunk(dtype: Dtype,
                           count: int = 1,
                           start: Union[None, int] = None) -> Union[Tuple[Any], Any]:
                if isinstance(start, int):
                    dat.seek(start)

                fmt = dtype2fmt[dtype]*count
                constructor = dtype2constr[dtype]
                dtype_size = dtype.value

                if dtype in (Dtype.wchar_t, Dtype.char):
                    return constructor(dat.read(dtype_size*count)).replace("\x00", "")
                else:
                    retn = struct.unpack("<" + fmt, dat.read(dtype_size * count))

                if len(retn) == 1:
                    return constructor(retn[0])
                else:
                    return tuple(constructor(i) for i in retn)

            cSignature = read_chunk(Dtype.char, 4)

            # Read the APT file header --------------------------------------------------------------------------------
            iHeaderSize = read_chunk(Dtype.int32)
            iHeaderVersion = read_chunk(Dtype.int32)
            wcFileName = read_chunk(Dtype.wchar_t, 256)
            ftCreationTime = read_chunk(Dtype.filetime)
            llIonCount = read_chunk(Dtype.int64)

            if verbose:
                print(f"\nReading header of {filepath}")
                print(f"\tcSignature: " + cSignature)
                print(f"\tiHeaderSize: {iHeaderSize}")
                print(f"\tiHeaderVersion: {iHeaderVersion}")
                print(f"\twcFileName: {wcFileName}")
                print(f"\tftCreationTime: {ftCreationTime}")
                print(f"\t11IonCount: {llIonCount}")

            # Read the APT sections ----------------------------------------------------------------------------
            section_start = iHeaderSize
            section_data = {}

            while True:
                sec_sig = read_chunk(Dtype.char, 4, section_start)

                if sec_sig == '':
                    # EOF reached
                    break

                # Flag used to not include a section in the Roi when a configuration
                # situation is not implemented or handled
                skip_sec = False

                sec_header_size = read_chunk(Dtype.int32)
                sec_header_ver = read_chunk(Dtype.int32)
                sec_type = read_chunk(Dtype.wchar_t, 32)
                sec_ver = read_chunk(Dtype.int32)

                sec_rel_type = RelType(read_chunk(Dtype.int32))
                is_one_to_one = sec_rel_type == RelType.ONE_TO_ONE
                if not is_one_to_one:
                    warn(f"APAV does not handle REL_TYPE != ONE_TO_ONE, section \"{sec_type}\" will be ignored")
                    skip_sec = True

                sec_rec_type = RecordType(read_chunk(Dtype.int32))
                is_fixed_size = sec_rec_type == RecordType.FIXED_SIZE
                if not is_fixed_size:
                    warn(f"APAV does not handle RECORD_TYPE != FIXED_SIZE, section \"{sec_type}\" will be ignored")
                    skip_sec = True

                sec_rec_dtype = RecordDataType(read_chunk(Dtype.int32))
                if sec_rec_dtype in (RecordDataType.DT_UNKNOWN, RecordDataType.OTHER, RecordDataType.CHARSTRING):
                    warn(f"APAV does not handle RECORD_TYPE == {sec_rec_dtype}, section \"{sec_type}\" will be ignored")
                    skip_sec = True

                sec_dtype_size = read_chunk(Dtype.int32)
                sec_rec_size = read_chunk(Dtype.int32)
                sec_data_unit = read_chunk(Dtype.wchar_t, 16)
                sec_rec_count = read_chunk(Dtype.int64)
                sec_byte_count = read_chunk(Dtype.int64)

                if verbose:
                    print("\nReading new section")
                    print(f"\tSection header sig: {sec_sig}")
                    print(f"\tSection header size: {sec_header_size}")
                    print(f"\tSection header version: {sec_header_ver}")
                    print(f"\tSection type: {sec_type}")
                    print(f"\tSection version: {sec_ver}")
                    print(f"\tSection relative type: {sec_rel_type}")
                    print(f"\tSection record type: {sec_rec_type}")
                    print(f"\tSection record data type: {sec_rec_dtype}")
                    print(f"\tSection data type size (bits): {sec_dtype_size}")
                    print(f"\tSection record size: {sec_rec_size}")
                    print(f"\tSection data type unit: {sec_data_unit}")
                    print(f"\tSection record count: {sec_rec_count}")
                    print(f"\tSection byte count: {sec_byte_count}")

                if not skip_sec:
                    columns = int(sec_rec_size/(sec_dtype_size/8))
                    records = int(sec_rec_count)
                    count = records*columns
                    in_data = np.fromfile(filepath,
                                         record_dtype2numpy_dtype(sec_rec_dtype, sec_dtype_size),
                                         count,
                                         offset=section_start+sec_header_size)
                    if columns > 1:
                        section_data[sec_type] = in_data.reshape(records, columns)
                    else:
                        section_data[sec_type] = in_data

                section_start = section_start + sec_byte_count + sec_header_size

        has_mass_data = "Mass" in section_data.keys()
        has_pos_data = "Position" in section_data.keys()

        # Map some APT section names to names that APAV expects, otherwise the provided name is retained
        name_map = {
            "Multiplicity": "ipp",
            "Time of Flight": "tof",
            "XDet_mm": "det_x",
            "YDet_mm": "det_y",
            "Voltage": "dc_voltage",
            "Pulse Voltage": "pulse_voltage"
        }

        # Require mass and position data, clean up some sections, and account for possible duplicate sections (i.e.
        # XDet_mm + YDet_mm combined with Detector Coordinates
        if not has_mass_data:
            raise AttributeError("APT file must have include a mass section")
        elif not has_pos_data:
            raise AttributeError("APT file must have include a position section")

        mass = section_data.pop("Mass")
        pos = section_data.pop("Position")
        data = np.concatenate((pos,mass.reshape(-1,1)),axis=1)

        # There are 2 difference ways that detector space coordinates can be included in an apt file, as a single
        # section containing both x/y or the x and y in separate sections. Only when the separate x/y sections are not
        # present we will load the combined x/y data (which we separate into different x and y arrays).
        if "Detector Coordinates" in section_data.keys():
            temp = section_data.pop("Detector Coordinates")
            if "XDet_mm" not in section_data.keys():
                section_data["det_x"] = temp[:, 0]
            if "YDet_mm" not in section_data.keys():
                section_data["det_y"] = temp[:, 1]

        #aptdata = cls(pos, mass)
        #aptdata._filepath = filepath
        #for key, data in section_data.items():
        #    name = key if key not in name_map.keys() else name_map[key]
        #    aptdata.misc[name] = data

        return pos,mass


#-------------------------------------------------------------------------------
# TEST
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    data_3D = np.array([[1,2,3, 0], [4,5,6, 0], [7,8,9, 0]])
    point_3D = AtomCloud()
    point_3D.init_point_cloud(data_3D[:,0:3],data_3D[:,3])
    dist, ind = point_3D.query_radius(0, 20)
    print(dist)
    print(ind)