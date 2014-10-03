"""
.. moduleauthor:: Mathieu Westphal < mathieu.westphal@obs.ujf-grenoble.fr >


"""
from libc.stdlib cimport malloc, free
from libc.string cimport strcmp
from cpython.string cimport PyString_AsString

cdef extern from "../../src/libastrochem.h":
    cdef double CHI_DEFAULT
    cdef double COSMIC_DEFAULT
    cdef double GRAIN_SIZE_DEFAULT
    cdef double ABS_ERR_DEFAULT
    cdef double REL_ERR_DEFAULT

    ctypedef struct params_t:
        double nh
        double tgas
        double tdust

    ctypedef struct phys_t:
        double chi
        double cosmic
        double grain_size
        double grain_abundance

    ctypedef struct net_t:
        int n_species
        char **species_names

    ctypedef struct cell_t:
        double av
        double nh
        double tgas
        double tdust

    ctypedef struct astrochem_mem_t:
        params_t params

    void read_network ( const char *chem_file, net_t* network, const int verbose)
    void free_network ( net_t* network)
    int alloc_abundances( const net_t* network, double** abundances )
    void free_abundances( double* abundances )
    int set_initial_abundances( char** species, int n_initialized_abundances, const double* initial_abundances, const net_t* network, double* abundances )
    int solver_init( const cell_t* cell, const net_t* network, const phys_t* phys, const double* abundances , double density, double abs_err, double rel_err, astrochem_mem_t* astrochem_mem )
    int solve( const astrochem_mem_t* astrochem_mem, const net_t* network, double* abundances, double time , const cell_t* new_cell, int verbose )
    void solver_close( astrochem_mem_t* astrochem_mem )

_ABS_ERR_DEFAULT = ABS_ERR_DEFAULT
_REL_ERR_DEFAULT = REL_ERR_DEFAULT

cdef class Network:
    """
    Chemical Network class

    :param chem_file:  Chemical network file string to load network from and use in solver.
    :type chem_file: const char*
    :param verbose: verbose if 1, quiet if 0
    :type verbose: int

    """
    cdef net_t thisstruct

    #for sphinx doc only
    def __init__( self, const char* chem_file, int verbose ):
        """
        __init__( chem_file, verbose )
        """

    #real cython init
    def __cinit__( self, const char* chem_file, int verbose ):
        read_network( chem_file, &self.thisstruct, verbose )

    def __dealloc__(self):
        free_network( &self.thisstruct )

cdef class Phys:
    """
    Physical parameters to use in chemical reaction solver
    """
    cdef public phys_t thisstruct

    #for Sphinx doc only
    def __init__( self ):
        """
        __init__( )
        """

    #real cython init
    def __cinit__( self ):
        self.thisstruct.chi = CHI_DEFAULT
        self.thisstruct.cosmic = COSMIC_DEFAULT
        self.thisstruct.grain_size = GRAIN_SIZE_DEFAULT
        self.thisstruct.grain_abundance = 0

    property chi:
        """
        Chi physical property

        :getter: Returns this chi physical property
        :setter: Sets this chi physical property
        :type: double
        """
        def __get__(self):
            return self.thisstruct.chi
        def __set__(self, double chi):
            self.thisstruct.chi = chi

    property cosmic:
        """
        Cosmic physical property

        :getter: Returns this cosmic physical property
        :setter: Sets this cosmic physical property
        :type: double
        """
        def __get__(self):
            return self.thisstruct.cosmic
        def __set__(self, double cosmic):
            self.thisstruct.cosmic = cosmic

    property grain_size:
        """
        Grain Size physical property

        :getter: Returns this grain size physical property
        :setter: Sets this grain size physical property
        :type: double
        """
        def __get__(self):
            return self.thisstruct.grain_size
        def __set__(self, double grain_size):
            self.thisstruct.grain_size = grain_size

    property grain_abundance:
        """
        Grain Abundance physical property

        :getter: Returns this grain abundance physical property
        :setter: Sets this grain abundance physical property
        :type: double
        """
        def __get__(self):
            return self.thisstruct.grain_abundance
        def __set__(self, double grain_abundance):
            self.thisstruct.grain_abundance = grain_abundance

cdef class Cell:
    """
    Chemical cell class

    :param av:  av parameters of the cell.
    :type av: double
    :param nh:  nh parameters of the cell.
    :type nh: double
    :param tgas:  tgas parameters of the cell.
    :type tgas: double
    :param tdust:  tdust parameters of the cell.
    :type tdust: double
    """
    cdef public cell_t thisstruct

    #for Sphinx doc only
    def __init__( self, double av, double nh, double tgas, double tdust ):
        """
        __init__( av, nh, tgas, tdust )
        """

    # real cython init
    def __cinit__( self, double av, double nh, double tgas, double tdust ):
        self.thisstruct.av = av
        self.thisstruct.nh = nh
        self.thisstruct.tdust = tgas
        self.thisstruct.tgas = tdust

    property av:
        """
        av cell property

        :getter: Returns this av cell property
        :setter: Sets this av cell property
        :type: double
        """
        def __get__(self):
            return self.thisstruct.av
        def __set__(self, double av):
            self.thisstruct.av = av

    property nh:
        """
        nh cell property

        :getter: Returns this nh cell property
        :setter: Sets this nh cell property
        :type: double
        """
        def __get__(self):
            return self.thisstruct.nh
        def __set__(self, double nh):
            self.thisstruct.nh = nh

    property tgas:
        """
        tgas cell property

        :getter: Returns this tgas cell property
        :setter: Sets this tgas cell property
        :type: double
        """
        def __get__(self):
            return self.thisstruct.tgas
        def __set__(self, double tgas):
            self.thisstruct.tgas = tgas

    property tdust:
        """
        tdust cell property

        :getter: Returns this tdust cell property
        :setter: Sets this tdust cell property
        :type: double
        """
        def __get__(self):
            return self.thisstruct.tdust
        def __set__(self, double tdust):
            self.thisstruct.tdust = tdust

cdef class Solver:
    """
    Chemical solver class. Compute abundances in a network at a certain time using inital abundances and differents parameters.

    :param cell:  Chemical cell class to use in solver.
    :type cell: :class:`.Cell`
    :param chem_file:  Chemical network file string to load network from and use in solver.
    :type chem_file: const char*
    :param phys:  Physical properties class to use in solver.
    :type phys: :class:`.Phys`
    :param abs_err: Absolute acceptable error to use in solver.
    :type abs_err: double
    :param rel_err: Relative acceptable error to use in solver.
    :param initial_abundances: Initial abundances dictionnary with format {Species:Value}
    :type initial_abundances: dictionnary
    :param density: Density to use in solver
    :type denisty: double
    :param verbose: verbose if 1, quiet if 0
    :type verbose: int
    :returns:  :class:`.Solver`
    :raises: MemoryError

    """

    cdef astrochem_mem_t astrochemstruct
    cdef double* abundances
    cdef Network network
    cdef int verbose

    # For Sphinx only (only used for docstring)
    def __init__( self , cell , chem_file, phys, abs_err, rel_err, initial_abundances , density, verbose ):
        """
        __init__( cell , chem_file, phys, abs_err, rel_err, initial_abundances , density, verbose )
        """

    # Real Cython init
    def __cinit__( self , cell , const char* chem_file, phys, abs_err, rel_err, initial_abundances , double density, int verbose ):
        self.verbose = verbose
        self.network = Network( chem_file, verbose )
        cdef net_t c_net = self.network.thisstruct
        cdef phys_t c_phys = phys.thisstruct
        cdef cell_t c_cell = cell.thisstruct

        if alloc_abundances( &c_net , &self.abundances ) != 0 :
            raise MemoryError
        cdef char **initial_abundances_str = <char **>malloc(len(initial_abundances) * sizeof(char *))
        cdef double* initial_abundances_val = <double*>malloc(len(initial_abundances) * sizeof(double))
        cdef int j = 0
        for i in initial_abundances:
            initial_abundances_str[j] = PyString_AsString(i)
            initial_abundances_val[j] = initial_abundances[i]
            j+=1
        set_initial_abundances( < const char** >initial_abundances_str, len(initial_abundances), initial_abundances_val,  &c_net, self.abundances )
        free( initial_abundances_str )
        free( initial_abundances_val )

        solver_init( &c_cell, &c_net, &c_phys , self.abundances, density, abs_err, rel_err, &self.astrochemstruct )

    def __dealloc__(self):
        free_abundances( self.abundances )
        solver_close(  &self.astrochemstruct )

    def solve(self, time , new_cell=0 ):
        """
        Solve chemical reaction for a certain time

        :param time:  time to solve the system at.
        :type time: double
        :param new_cell:  cell class to use in solver, optionnal
        :type new_cell: :class:`.Cell`
        :returns:  dictionnary -- computed abundances.
        :raises: ArithmeticError
        """
        cdef net_t c_net = self.network.thisstruct
        cdef cell_t c_new_cell
        if new_cell== 0 :
            if solve( &self.astrochemstruct, &c_net , self.abundances, time , NULL, self.verbose ) != 0:
                raise ArithmeticError("Solve Error")
        else:
            c_new_cell = new_cell.thisstruct
            if solve( &self.astrochemstruct, &c_net , self.abundances, time , &c_new_cell, self.verbose ) != 0:
                raise ArithmeticError("Solve Error")
        cdef int i
        cdef bytes py_string
        ret = {}
        for i in range( c_net.n_species ):
            py_string =  c_net.species_names[i]
            ret[py_string] = self.abundances[i]
        return ret
