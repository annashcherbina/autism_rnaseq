# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );

# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0050886","endocrine process",0.276196381263739,-0.860943681319978,1.50289854629499,1.69897000433602,1E-300,0.934875028333895,0),
                     c("GO:1902236","negative regulation of endoplasmic reticulum stress-induced intrinsic apoptotic signaling pathway",0.112733216842343,6.43978865394254,2.49794610903861,1.32221929473392,1E-300,0.878829968359696,0),
                     c("GO:1903513","endoplasmic reticulum to cytosol transport",0.0958232343159912,-0.29888654492626,-6.32371441395609,1.25527250510331,1E-300,0.943401775208963,0),
                     c("GO:0071888","macrophage apoptotic process",0.0169099825263514,1.76964688945208,6.99616235570034,0.602059991327962,1E-300,0.997167302123937,0.00297791),
                     c("GO:0051301","cell division",2.78451045600586,3.2776435718254,6.32146753898226,2.69460519893357,1E-300,0.995787749252727,0.0044211),
                     c("GO:0098743","cell aggregation",0.11836987768446,-1.75001651826103,6.55376190128351,1.34242268082221,0.000477991850914072,0.996758569876336,0.00451162),
                     c("GO:0051298","centrosome duplication",0.202919790316217,2.47336341575932,-7.54084038132231,1.56820172406699,1E-300,0.895029780642427,0.00476024),
                     c("GO:0042403","thyroid hormone metabolic process",0.11836987768446,-0.33506763096806,-2.66420295028861,1.34242268082221,1E-300,0.928398468980089,0.01639409),
                     c("GO:1900006","positive regulation of dendrite development",0.0958232343159912,5.49274579645006,-4.2528149717794,1.25527250510331,1E-300,0.897663765210266,0.02520546),
                     c("GO:0032465","regulation of cytokinesis",0.512936136632659,-3.36996205693621,6.33973816761682,1.96378782734556,1E-300,0.975724119886796,0.02871156),
                     c("GO:1902692","regulation of neuroblast proliferation",0.186009807789865,0.338065049674761,5.16208712476281,1.53147891704226,0.000876685514579726,0.951002443993644,0.02995534),
                     c("GO:0051983","regulation of chromosome segregation",0.541119440843244,4.78839539105002,5.69505383959774,1.98677173426624,0.000860817903456111,0.975615834625902,0.03300408),
                     c("GO:0042362","fat-soluble vitamin biosynthetic process",0.0394566258948199,-6.44059651856366,-3.62755348896238,0.903089986991944,1E-300,0.952036605251553,0.075682),
                     c("GO:0006740","NADPH regeneration",0.0845499126317569,-0.139734949644504,7.25493173348391,1.20411998265592,1E-300,0.976652793333226,0.07985728),
                     c("GO:0045821","positive regulation of glycolytic process",0.112733216842343,2.8719763610478,3.5600690536936,1.32221929473392,1E-300,0.958294471083811,0.11474372),
                     c("GO:0035637","multicellular organismal signaling",0.665125979369821,-3.52762164788,3.94621523710055,2.07554696139253,1E-300,0.947778700940774,0.11502594),
                     c("GO:0060973","cell migration involved in heart development",0.0958232343159912,-5.12042571161228,-0.776078725842117,1.25527250510331,1E-300,0.88306072023609,0.12659773),
                     c("GO:0060457","negative regulation of digestive system process",0.0901865734738741,5.29949874070603,-0.786594912625077,1.23044892137827,1E-300,0.917285678944747,0.12660131),
                     c("GO:0035970","peptidyl-threonine dephosphorylation",0.0901865734738741,-4.05336000895381,-4.92396033195835,1.23044892137827,1E-300,0.959352739159886,0.1360348),
                     c("GO:0070493","thrombin-activated receptor signaling pathway",0.0507299475790542,-4.87181260118756,5.11276967886932,1,0.000710255402511067,0.970183654103693,0.14234157),
                     c("GO:0007566","embryo implantation",0.219829772842568,-6.08270152249893,2.48576037418447,1.60205999132796,1E-300,0.939734296640844,0.15722417),
                     c("GO:0099587","inorganic ion import across plasma membrane",0.507299475790542,-1.55327455046858,-6.40429963704687,1.95904139232109,1E-300,0.935591029472266,0.19789345),
                     c("GO:0021988","olfactory lobe development",0.1972831294741,-5.51308703806612,1.37310570960086,1.55630250076729,1E-300,0.900948801914996,0.21936919),
                     c("GO:0030206","chondroitin sulfate biosynthetic process",0.101459895158108,-5.22732835939822,-4.4660418801379,1.27875360095283,1E-300,0.919642222855195,0.24215005),
                     c("GO:0071459","protein localization to chromosome, centromeric region",0.11836987768446,0.0768822070837475,-7.05812010850824,1.34242268082221,1E-300,0.952630514237082,0.25670473),
                     c("GO:0001672","regulation of chromatin assembly or disassembly",0.180373146947748,4.63702010003945,-5.87360488736628,1.51851393987789,1E-300,0.922297129981222,0.27019495),
                     c("GO:0061035","regulation of cartilage development",0.383292937263965,6.24425322569884,-2.58855459729641,1.83884909073726,1E-300,0.920773733763445,0.27865119),
                     c("GO:0042310","vasoconstriction",0.180373146947748,0.176980022782886,-0.284963765447425,1.51851393987789,1E-300,0.896559968083094,0.31219605),
                     c("GO:0045665","negative regulation of neuron differentiation",0.405839580632433,7.10975687724962,-1.39917760208092,1.86332286012046,1E-300,0.90612873082239,0.31498183),
                     c("GO:0007588","excretion",0.248013077053154,-0.306653714000672,1.78924353348444,1.65321251377534,1E-300,0.935338353423087,0.32055681),
                     c("GO:1901028","regulation of mitochondrial outer membrane permeabilization involved in apoptotic signaling pathway",0.11836987768446,1.68147375910061,-5.44033784765247,1.34242268082221,1E-300,0.809662003868349,0.34755626),
                     c("GO:0040037","negative regulation of fibroblast growth factor receptor signaling pathway",0.101459895158108,6.85121615145673,2.09963018055355,1.27875360095283,1E-300,0.903358755278654,0.36685604),
                     c("GO:0050777","negative regulation of immune response",0.986415647370498,6.68524264472506,1.31103476839986,2.24551266781415,1E-300,0.862602199417409,0.38929646),
                     c("GO:0006002","fructose 6-phosphate metabolic process",0.0620032692632884,-4.83948980094723,-5.33136000419259,1.07918124604762,0.00087266402139358,0.948891624678205,0.3941751),
                     c("GO:0045651","positive regulation of macrophage differentiation",0.0845499126317569,6.15565345436557,-1.78202903892455,1.20411998265592,1E-300,0.875706543273547,0.39797729),
                     c("GO:0060055","angiogenesis involved in wound healing",0.0958232343159912,-6.40766054486624,-0.017390782124048,1.25527250510331,1E-300,0.91228785360771,0.41144675),
                     c("GO:0061041","regulation of wound healing",0.710219266106758,6.25001553334145,3.3592488872964,2.10380372095596,1E-300,0.903996531008509,0.41262763),
                     c("GO:0070839","metal ion export",0.236739755368919,-1.68920673344615,-6.84906908049103,1.63346845557959,1E-300,0.940598616685912,0.44348473),
                     c("GO:0048557","embryonic digestive tract morphogenesis",0.0958232343159912,-6.51765987663857,0.608427708964563,1.25527250510331,1E-300,0.916384997182357,0.45736176),
                     c("GO:0033044","regulation of chromosome organization",1.18933543768671,4.11963254680781,-5.73996728415418,2.32633586092875,1E-300,0.926585709572959,0.45932768),
                     c("GO:0003214","cardiac left ventricle morphogenesis",0.0845499126317569,-5.97314443607528,0.16432061874986,1.20411998265592,1E-300,0.909463482910088,0.48507415),
                     c("GO:0010460","positive regulation of heart rate",0.135279860210811,3.48930806917739,-1.80457904943489,1.39794000867204,1E-300,0.904647962296536,0.49584321));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$log_size <- as.numeric( as.character(one.data$log_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);


# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("/path_to_your_file/revigo-plot.pdf");

