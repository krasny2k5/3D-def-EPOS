Data for PAF inversion for the Lorca Area.

There are two folders:

- Originals: The original CPT processing with an extended area.
- Originals - ROI: The data from the previous folder with a reduced extend (used in the paper)

Inside each folder there is a .sh program parser_LOS(_original).sh which converts the data to the PAF format. It uses AWK and it is intended to be run in a linux machine.

The coverted data is inside "asc" and "dsc" folders of each type (full extend and ROI clipped). The dates used are the following:

ASCENDING (22 dates)
12/11/2015 06/12/2015 30/12/2015 23/01/2016 16/02/2016 11/03/2016 04/04/2016 28/04/2016	10/05/2016 03/06/2016 27/06/2016 21/07/2016 14/08/2016 07/09/2016 01/10/2016 25/10/2016 06/11/2016 30/11/2016	24/12/2016 05/01/2017 29/01/2017 22/02/2017


DESCENDING (19 dates)
12/12/2015 24/12/2015 17/01/2016 10/02/2016 05/03/2016 29/03/2016 22/04/2016 16/05/2016 09/06/2016 03/07/2016 27/07/2016 20/08/2016 13/09/2016 07/10/2016 31/10/2016 24/11/2016 18/12/2016 11/01/2017	28/02/2017

The files are numbered as D0000**.txt where ** is the corresponding number using an ascending date criteria (D000001.txt is the first date, D000002.txt second, etc...)

The data is distributed across 9 different columns:

UTM_EASTING UTM_NORTHING HEIGHT LOS_VALUE 0 0 1 0 0 for example:

625357.38 4175143.75 316.29 0 0 0 1 0 0
625578.44 4175065.5 313.29 0 0 0 1 0 0
625567.12 4175025 316.89 0 0 0 1 0 0
624939.44 4175084.5 320.74 0 0 0 1 0 0
625388.31 4174969.25 314.63 0 0 0 1 0 0
624809.88 4175063.25 329.32 0 0 0 1 0 0
625769.06 4174865 318.34 0 0 0 1 0 0


Last 6 columns correspond with displacement (dz, dx, dy) and the weight for each value (between 0 and 1, 0 must be used for bad data and 1 for good data). The decimal separator must be dot "." and the column separator must be a space in order to be processed by PAF software.

This dataset is the same used in the 2018 Scientific Reports article: Modeling the two- and three-dimensional displacement field in Lorca, Spain, subsidence and the global implications (https://doi.org/10.1038/s41598-018-33128-0)
