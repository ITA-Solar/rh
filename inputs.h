/* ------- file: -------------------------- inputs.h ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Jul 24 12:38:33 2013 --

       --------------------------                      ----------RH-- */

#ifndef __INPUTS_H__
#define __INPUTS_H__


/* --- Include file for command line options and keyword input
       parameters. --                                  -------------- */


#define MAX_KEYWORD_LENGTH  32
#define MAX_VALUE_LENGTH    80

#define CHECK_OPTION(option, name, length) \
        (strncmp(option, name, MAX(((int) length), strlen(option))) == 0)


enum keywordtype  {KEYWORD_REQUIRED, KEYWORD_DEFAULT, KEYWORD_OPTIONAL};
enum S_interpol   {S_LINEAR, CUBIC_HERMITE, BEZIER, BEZIER3};
enum S_interpol_stokes   {DELO_PARABOLIC, DELO_BEZIER3};
enum order_3D     {LINEAR_3D, BICUBIC_3D};
enum ne_solution  {NONE, ONCE, ITERATION};


typedef struct {
  char   name[MAX_KEYWORD_LENGTH];
  int    minlength;
  bool_t value_required;
  char   value[MAX_VALUE_LENGTH];
  void  *pointer, (*setValue)(char *value, void *pointer);
  char   message[MAX_VALUE_LENGTH];
} Option;

typedef struct {
  char   keyword_input[MAX_VALUE_LENGTH], *wavetable, *molecule;
  bool_t quiet, showkeywords;
  FILE  *logfile;
} CommandLine;

typedef struct {
  char   keyword[MAX_KEYWORD_LENGTH], value[MAX_VALUE_LENGTH];
  bool_t set;
  enum   keywordtype type;
  void  *pointer, (*setValue)(char *value, void *pointer);
} Keyword;

typedef struct {
  char   atmos_input[MAX_VALUE_LENGTH],
         abund_input[MAX_VALUE_LENGTH],
         wavetable_input[MAX_VALUE_LENGTH],
         atoms_input[MAX_VALUE_LENGTH],
         molecules_input[MAX_VALUE_LENGTH],
         Stokes_input[MAX_VALUE_LENGTH],
         KuruczData[MAX_VALUE_LENGTH],
         pfData[MAX_VALUE_LENGTH],
         BarklemDir[MAX_VALUE_LENGTH],
         fudgeData[MAX_VALUE_LENGTH],
         atmos_output[MAX_VALUE_LENGTH],
         spectrum_output[MAX_VALUE_LENGTH],
         geometry_output[MAX_VALUE_LENGTH],
         opac_output[MAX_VALUE_LENGTH],
         JFile[MAX_VALUE_LENGTH],
         background_File[MAX_VALUE_LENGTH],
         background_ray_File[MAX_VALUE_LENGTH],
         H_atom[MAX_VALUE_LENGTH],
         H2_molecule[MAX_VALUE_LENGTH],
         radrateFile[MAX_VALUE_LENGTH],
         collrateFile[MAX_VALUE_LENGTH],
         dampingFile[MAX_VALUE_LENGTH],
         coolingFile[MAX_VALUE_LENGTH],
         Itop[MAX_VALUE_LENGTH];
  bool_t magneto_optical, XRD, Eddington,
         backgr_pol, limit_memory, allow_passive_bb, NonICE,
         rlkscatter, prdh_limit_mem, backgr_in_mem, xdr_endian,
         old_background, accelerate_mols;
  enum   solution startJ;
  enum   StokesMode StokesMode;
  enum   S_interpol S_interpolation;
  enum   S_interpol_stokes S_interpolation_stokes;
  enum   order_3D interpolate_3D;
  enum   ne_solution solve_ne;
  enum   PRDangle PRD_angle_dep;
  int    isum, Ngdelay, Ngorder, Ngperiod, NmaxIter,
         PRD_NmaxIter, PRD_Ngdelay, PRD_Ngorder, PRD_Ngperiod,
         NmaxScatter, Nthreads;
  /* Tiago, for collisional-radiative switching */
  double crsw, crsw_ini;
  /* Tiago, for PRD switching */
  double prdswitch, prdsw;
  /* Tiago, for micro turbulence multiplication and addition */
  double vturb_mult, vturb_add;
  /* Tiago, for escape probability iterations */
  int    NpescIter;
  /* Tiago, added this for 1.5D version */
  int    p15d_nt, p15d_x0, p15d_x1, p15d_xst, p15d_y0, p15d_y1, p15d_yst;
  int    Natoms;
  double p15d_tmax;
  bool_t p15d_wxtra, p15d_rerun, p15d_refine, p15d_zcut, p15d_wtau;
  bool_t p15d_wpop, p15d_wrates, p15d_wcrates;
  int    p15d_flush_interval;
  double iterLimit, PRDiterLimit, metallicity, *wavetable;
  unsigned int Nxwave;
  /* Tiago, for saving the input files */
  char  *atoms_file_contents, *keyword_file_contents, *ray_file_contents;
  char **atomic_file_contents;
  char  *kurucz_file_contents, **kurucz_line_file_contents;
  char **kurucz_line_file_name;
  int Nkurucz_files;
  pthread_attr_t thread_attr;
  /* --- Node-partitioned pool mode (node-shared atmosphere cache) ----
     Added at the END of the struct so existing field offsets are
     preserved across incremental builds — old object files using a
     stale inputs.h still find atoms_file_contents etc. at the right
     place.  When TRUE (default), rh15d_ray_pool reads the atmosphere
     once per node into a shared-memory window and serves every
     readAtmos call from that cache instead of touching the file.
     atmos_cache_cyclic: v1 uses contiguous row-based decomposition.
     When TRUE, use cyclic row decomposition for spatial load balance.
     Plumbed but not yet honored in v1. */
  bool_t use_node_atmos_cache;
  bool_t atmos_cache_cyclic;
  /* Pool mode: when TRUE (default), writeCollective_pool uses
     H5FD_MPIO_COLLECTIVE for the phase-2 flush; when FALSE, each rank
     writes its own hyperslabs with H5FD_MPIO_INDEPENDENT.  Independent
     mode is useful on filesystems where collective buffering adds more
     overhead than it saves, or for debugging MPI-IO behaviour. */
  bool_t pool_collective_write;
  /* Lustre contention-reduction toggle: when TRUE (default), the
     create_hdf5_fapl helpers pass the romio/stripe/cb_nodes MPI-IO
     hints (built in initParallel) to H5Pset_fapl_mpio, AND enable
     collective HDF5 metadata ops on the output FAPL.  When FALSE,
     MPI_INFO_NULL is passed and collective metadata ops are left off
     — toggle off to A/B against pre-42fc036/e1767d1 baseline.
     Note: H5Pset_alignment and H5Pset_meta_block_size are always
     applied (pure HDF5-level tuning with no behavioural side effects);
     H5Pset_file_locking(false) is also always applied (correctness
     requirement on Lustre, not a performance knob). */
  bool_t mpiio_lustre_hints;
} InputData;


/* --- Associated function prototypes --               -------------- */

int   getLine(FILE *inputFile, char *commentChar, char *line,
	      bool_t exit_on_EOF);
int   getLineString(char **inputString, char *commentChar, char *line,
                    bool_t exit_on_EOF);
void  parse(int argc, char *argv[], int Noption, Option *theOptions);
void  readInput(char *input_string);
void  readValues(char *fp_keyword, int Nkeyword, Keyword *theKeywords);

char *readWholeFile(const char *filename);
char *sgets(char *str, int num, char **input);

void  setAngleSet(char *value, void *pointer);
void  setcharValue(char *value, void *pointer);
void  setboolValue(char *value, void *pointer);
void  setintValue(char *value, void *pointer);
void  setdoubleValue(char *value, void *pointer);
void  setstartValue(char *value, void *pointer);
void  setnesolution(char *value, void *pointer);
void  setStokesMode(char *value, void *pointer);
void  setPRDangle(char *value, void *pointer);
void  setThreadValue(char *value, void *pointer);
void  setInterpolate_3D(char *value, void *pointer);
void  set_S_interpolation(char *value, void *pointer);
void  set_S_interpolation_stokes(char *value, void *pointer);
void  showValues(int Nkeyword, Keyword *theKeywords);
char *substring(const char *string, int N0, int Nchar);
void  UpperCase(char *string);

void  writeInput();

#endif /* !__INPUTS_H__ */

/* ------- end ---------------------------- inputs.h ---------------- */
