/**
* Name: tissuedetectorv2
* Author: Sandra Saez
* Description: 
*/

model tissuedetectorv2

global torus:false {
	// Initial number of robots
	int nb_robots_init;
	int grid_size <- 50;
	tissue_cell my_cell_ini;
	
	// Diffusion  rate from input
	float diff_rate_chem1;
	float diff_rate_chem2;
	float diff_rate_chem3;
	
	// Emision rate of chemicals.
	float generation_rt_chem1;
	float generation_rt_chem2;
	float generation_rt_chem3;
	
	// Chemical 1 threshold
	float chem1_threshold <- 0.1;
	
	// Probabilities od differentiation
	float prob_differentiation <- 0.5;
	
	// Sensing threshold
	float chem1_sensing_th <- 0.0000001;
	float chem2_sensing_th <- 0.0000001;
	float chem3_sensing_th <- 0.0000001;
	
	// Evaporation of chemicals
	float evaporation_per_cycle <- 0.01 min: 0.0 max: 2.0;
	
	// Tumor cell: division probability
	float div_prob <- 0.01;
	
	// Read image of the tissue
	file map_init <- image_file("../images/damaged_tissue3.png");
	
	// Set up the world
	init {
		matrix init_data <- map_init as_matrix {grid_size,grid_size};
		int j <- 0; // row
		int i <- 0; // column
		loop times:nb_robots_init+1 {
			if i < 10 {
				my_cell_ini <- tissue_cell[i+5,j+5];
				create robot with:[my_cell::my_cell_ini, location::my_cell_ini.location];
				i <- i + 1;
			} else {
				i <- 0;
				j <- j + 1;
			}
		}
		
		
		//create robot number:nb_robots_init;
		ask tissue_cell {
			loop cell over:init_data {
				rgb cell_color <- rgb(init_data[grid_x,grid_y]);
				if cell_color != #white {
					self.is_damaged <- true;
				}
			}
			
			// Create the damaged_cell agents
			if self.is_damaged {
				create damaged_cell number:1 with:[location::self.location, my_cell::tissue_cell(self.location)];
			}
		}	
	}
	
	reflex diffuse {
      diffuse var:chem1 on:tissue_cell proportion: diff_rate_chem1 radius:10 propagation: gradient;
      diffuse var:chem2 on:tissue_cell proportion: diff_rate_chem2 radius:10 propagation: gradient;
      diffuse var:chem3 on:tissue_cell proportion: diff_rate_chem3 radius:10 propagation: gradient;
    }
}

// Agent that represents all the explored tissue (a sample). The chemicals will be segregated into its cells.
grid tissue_cell height:grid_size width:grid_size neighbors:8 {
	bool is_damaged <- false;
	rgb color <- rgb(int(255 * (1 - chem1)), 255 * (1 - chem2), int(255 * (1 - chem3))) update:rgb(int(255 * (1 - chem1)), 255 * (1 - chem2), int(255 *(1 - chem3))) ;
	
	// Concentrations of the three different chemicals.
	float chem1 <- 0.0 update: (chem1<=evaporation_per_cycle) ? 0.0 : chem1-evaporation_per_cycle; 
	float chem2 <- 0.0 update: (chem2<=evaporation_per_cycle) ? 0.0 : chem2-evaporation_per_cycle;
	float chem3 <- 0.0 update: (chem3<=evaporation_per_cycle) ? 0.0 : chem3-evaporation_per_cycle;
	
	// List of neighbours at distance 1
	list<tissue_cell> neighbors <- self neighbors_at 1;
}

// Agent that represents each of the nanorobots.
species robot {
	float size <- 0.5;
	rgb color <- #lime;
	string chemicals <- "chem1" among: ["chem1", "chem2", "chem3"];
	string emitting <- "NO" among: ["NO", "CHEM1", "CHEM2", "CHEM3"];
	bool differenciation_checked <- false;
	
	tissue_cell my_cell update: tissue_cell(self.location);
	
	/* AUXILIARY FUNCTIONS */

	bool check_chem1 {
		return my_cell.chem1 > chem1_threshold;
	}
	
	// Check if there is some chemical in the current tissue_cell
	bool check_chemicals {
		return (my_cell.chem1 > chem1_sensing_th) and (my_cell.chem2 > chem2_sensing_th) and (my_cell.chem3 > chem1_sensing_th);
	}
	
	// If empty -> there is no tumor. Else -> there is tumor.
	bool check_tumor {
		return !empty(damaged_cell inside my_cell);
	}
	
	reflex emit when:(emitting != "NO"){
		if emitting = "CHEM1" {
			tissue_cell(my_cell).chem1 <- tissue_cell(my_cell).chem1+generation_rt_chem1;
		} else if emitting = "CHEM2" {
			tissue_cell(my_cell).chem2 <- tissue_cell(my_cell).chem2+generation_rt_chem2;
		} else if emitting = "CHEM3" {
			tissue_cell(my_cell).chem3 <- tissue_cell(my_cell).chem3+generation_rt_chem3;
		}
	}
	
	/*************************************** */
	/*             LIST OF ACTIONS                */
	/*************************************** */
	
	action random_walk {
		list<tissue_cell> my_cells_tmp <- my_cell.neighbors where (empty(robot inside each));
		if !empty(my_cells_tmp) {
			location <- one_of (shuffle(my_cells_tmp)).location;
		}	
	}
	
	action differentiate (string chemical){
		self.differenciation_checked <- true;
		if chemical = "chem1" {
			emitting <- "CHEM1";
		} else if chemical = "chem2" {
			emitting <- "CHEM2";
		} else if chemical = "chem3" {
			emitting <- "CHEM3";
		}
	}
	
	action kill_damaged_cell {
		damaged_cell choosen_cell <- one_of (damaged_cell inside my_cell);
		ask choosen_cell {
			self.my_cell.is_damaged <- false;
			do die;
		}
		do differentiate("chem1");
	}
	
	action move_up_gradient (string chemical) {
		list<tissue_cell> my_cells_tmp <- my_cell.neighbors where (empty(robot inside each));
		if !empty(my_cells_tmp) {
			if chemical = "chem1" {
				tissue_cell max_conc_cell <- my_cells_tmp with_max_of (each.chem1);
				location <- max_conc_cell.location;
			} else if chemical = "chem2" {
				tissue_cell max_conc_cell <- my_cells_tmp with_max_of (each.chem2);
				location <- max_conc_cell.location;
			} else if chemical = "chem3" {
				tissue_cell max_conc_cell <- my_cells_tmp with_max_of (each.chem3);
				location <- max_conc_cell.location;
			}
		}
	}
	
	/*************************************** */
	/*      DISPATCHER - BEHAVIOUR        */
	/*************************************** */
	
	
	reflex dispatcher when:!differenciation_checked {
		if (check_chemicals() = false) and (check_tumor() = false) {
			do random_walk();
		} else if (check_tumor() = true) {
			do kill_damaged_cell();
		} else if (my_cell.chem1 > chem1_threshold) {
			do random_walk();
		} else if (my_cell.chem1 > 0.0) {
			bool differentiated <- flip(prob_differentiation);
			if differentiated = true {
				do differentiate("chem2");
			} else {
				do move_up_gradient("chem1");
			}
		} else if (my_cell.chem2 > 0.0) {
			bool differentiated <- flip(prob_differentiation);
			if differentiated = true {
				do differentiate("chem3");
			} else {
				do move_up_gradient("chem2");
			}
		} else if (my_cell.chem3 > 0.0) {
			do move_up_gradient("chem2");
		}
	}
	
	aspect base {
		draw circle(size) color:color;
	}
}
	

species damaged_cell {
	float size <- 0.5;
	rgb color <- #red;
	tissue_cell my_cell;
	//tissue_cell my_cell <- tissue_cell(self.location);
	list<tissue_cell> recheable_cells;
	file my_icon <- file("../images/arnold_tuma.png");
	
	init {
		recheable_cells <- my_cell.neighbors where (empty(damaged_cell inside each));
	}
	
	/*reflex divide when:(recheable_cells != nil) and (flip(div_prob)) {
		tissue_cell my_cell_tmp <- one_of(recheable_cells);
		ask my_cell_tmp {
			create damaged_cell number:1 with:[location::self.location, my_cell::tissue_cell(self.location)];
		}
		
	}*/
	
	aspect base {
		draw circle(size) color:color;
	}
	
	aspect icon {
		draw my_icon size: 2*size ;
	}
}

experiment tissue_detector type: gui {
	parameter "Initial number of nanorobots: " var:nb_robots_init init:100 min:1 category:"Nanorobots";
	parameter "Emision rate of chemical 1" var:generation_rt_chem1 init:0.1 min: 0.0 category:"Nanorobots";
	parameter "Emision rate of chemical 2" var:generation_rt_chem2 init:0.1 min: 0.0 category:"Nanorobots";
	parameter "Emision rate of chemical 3" var:generation_rt_chem3 init:0.1 min: 0.0 category:"Nanorobots";
	parameter "Evaporation per cycle" var:evaporation_per_cycle init:0.01 min:0.0 category:"Chemicals";
	parameter "Diffusion rate of chemical 1" var:diff_rate_chem1 init:0.3 min: 0.0 max: 1.0 category:"Chemicals";
	parameter "Diffusion rate of chemical 2" var:diff_rate_chem2 init:0.3 min: 0.0 max: 1.0 category:"Chemicals";
	parameter "Diffusion rate of chemical 3" var:diff_rate_chem3 init:0.3 min: 0.0 max: 1.0 category:"Chemicals";
	parameter "Chemical 1 threshold" var:chem1_threshold init:0.2 min:0.0 category:"Chemicals";
	
	output {
		display main_display {
			grid tissue_cell lines:rgb("grey");
			species robot aspect:base;
			species damaged_cell aspect:base;
		}
	}
}