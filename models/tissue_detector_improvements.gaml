/**
* Name: tissuedetectorv2
* Author: Sandra Saez
* Description: 
*/

model tissuedetectorv2

global torus:false {
	//user_command create_robots_here action:create_robots with: [cell::tissue_cell(location)];
	
	// Inicialization method
	bool injected <- false;
	
	// Robots global variables
	int nb_robots_init;
	tissue_cell my_cell_ini;
	float max_t_differenciation;
	
	// Grid global variables
	int grid_size <- 40;
	
	// Damaged cell global variables
	float energy_reproduce;
	float nutrient_uptake;
	int nb_damaged_cells -> {length (damaged_cell)};
	
	// Diffusion  rate from input
	float diff_rate_chem1;
	float diff_rate_chem2;
	float diff_rate_chem3;
	
	// Emision rate of chemicals.
	float generation_rt_chem1;
	float generation_rt_chem2;
	float generation_rt_chem3;
	
	// Chemical 1 threshold
	float chem1_threshold <- 0.2;
	
	// Probabilities od differentiation
	float prob_differentiation <- 0.05;
	
	// Sensing threshold
	float chem1_sensing_th <- 0.0001;
	float chem2_sensing_th <- 0.0001;
	float chem3_sensing_th <- 0.0001;
	
	// Evaporation of chemicals
	float evaporation_per_cycle_chem1 <- 0.01 min: 0.0;
	float evaporation_per_cycle_chem2 <- 0.01 min: 0.0;
	float evaporation_per_cycle_chem3 <- 0.01 min: 0.0;
	
	// Tumor cell: division probability
	float div_prob <- 0.0;
	
	// Read image of the tissue
	file map_init <- image_file("../images/damaged_tissue3.png");
	
	// Set up the world
	init {
		matrix init_data <- map_init as_matrix {grid_size,grid_size};
		
		if injected = false {
			int j <- 0; // row
			int i <- 0; // column
			int side_dimension <- int(sqrt(world.nb_robots_init));
			loop times:nb_robots_init+1 {
				if i < side_dimension {
					my_cell_ini <- tissue_cell[i+10,j+10];
					create robot with:[my_cell::my_cell_ini, location::my_cell_ini.location, my_first_cell::tissue_cell(location)];
					i <- i + 1;
				} else {
					i <- 0;
					j <- j + 1;
				}
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
	float chem1 <- 0.0 update: (chem1<=evaporation_per_cycle_chem1) ? 0.0 : chem1-evaporation_per_cycle_chem1; 
	float chem2 <- 0.0 update: (chem2<=evaporation_per_cycle_chem2) ? 0.0 : chem2-evaporation_per_cycle_chem2;
	float chem3 <- 0.0 update: (chem3<=evaporation_per_cycle_chem3) ? 0.0 : chem3-evaporation_per_cycle_chem3;
	
	// List of neighbours at distance 1
	list<tissue_cell> neighbors <- self neighbors_at 1;
	
	action create_robots {
		create robot number:nb_robots_init with: [location::self.location, my_cell::tissue_cell(location), my_first_cell::tissue_cell(location)];
	}
	
	/*action create_robots_from_cell {
		int i <- 0;
		int j <- 0;
		int side_dimension <- int(sqrt(world.nb_robots_init));
		int init_i <- int(location.x);
		int init_j <- int(location.y);
		int original_column <- init_i;
		tissue_cell cell_tmp <- tissue_cell[init_i+i, init_j+j];
		loop times:nb_robots_init+1 {
			if i < side_dimension {	
				create robot with:[location::self.location, my_cell::cell_tmp];
				cell_tmp <- tissue_cell[init_i+i+1, init_i+j+1];
				i <- i+1;
			} else {
				j <- j+1;
				i <- original_column;
			}
		}
	} */
	
	user_command create_robots_here action:create_robots;
	//user_command create_robots_from_here action:create_robots_from_cell;
}

// Agent that represents each of the nanorobots.
species robot skills: [moving] {
	float size <- 0.5;
	rgb color <- #lime;
	string chemicals <- "chem1" among: ["chem1", "chem2", "chem3"];
	string emitting <- "NO" among: ["NO", "CHEM1", "CHEM2", "CHEM3"];
	string state <- "wander" among: ["wander", "follower", "emitter"];
	bool differenciation_checked <- false;
	
	// Temporal differenciation
	float max_time <- max_t_differenciation;
	float local_timer <- 0.0 max: max_time;
	
	tissue_cell my_cell update: tissue_cell(self.location);
	tissue_cell my_first_cell;
	
	/* AUXILIARY FUNCTIONS */

	bool check_chem1 {
		return my_cell.chem1 > chem1_threshold;
	}
	
	// Check if there is some chemical in the current tissue_cell
	bool check_chemicals {
		return (my_cell.chem1 > chem1_sensing_th or my_cell.chem2 > chem2_sensing_th or my_cell.chem3 > chem1_sensing_th);
	}
	
	// If empty -> there is no tumor. Else -> there is tumor.
	bool check_tumor {
		return !empty(damaged_cell inside my_cell);
	}
	
	reflex emit when:(emitting != "NO") {
		if emitting = "CHEM1" {
			tissue_cell(my_cell).chem1 <- tissue_cell(my_cell).chem1+generation_rt_chem1;
		} else if emitting = "CHEM2" {
			tissue_cell(my_cell).chem2 <- tissue_cell(my_cell).chem2+generation_rt_chem2;
		} else if emitting = "CHEM3" {
			tissue_cell(my_cell).chem3 <- tissue_cell(my_cell).chem3+generation_rt_chem3;
		}
		self.local_timer <- self.local_timer + 1.0;
	}
	
	reflex leave_emitting when:(local_timer >= max_time) and (state = "emitter") {
		self.state <- "wander";
		self.emitting <- "NO";
		self.local_timer <- 0.0; // Reset the timer
		self.differenciation_checked <- false; // To have the opportunity to differenciate again
	}
	
	reflex go_back_injection_site {
		if nb_damaged_cells = 0 {
			emitting <- "NO";
			do goto target:my_first_cell;  
		}
	}
	
	/*************************************** */
	/*             LIST OF ACTIONS                */
	/*************************************** */
	
	action random_walk {
		state <- "wander";
		list<tissue_cell> my_cells_tmp <- my_cell.neighbors where (empty(robot inside each));
		if !empty(my_cells_tmp) {
			location <- one_of (shuffle(my_cells_tmp)).location;
		}	
	}
	
	action differentiate (string chemical) {
		state <- "emitter";
		if chemical = "chem1" {
			emitting <- "CHEM1";
		} else if chemical = "chem2" {
			emitting <- "CHEM2";
		} else if chemical = "chem3" {
			emitting <- "CHEM3";
		}
	}
	
	action kill_damaged_cell {
		state <- "emitter";
		damaged_cell choosen_cell <- one_of (damaged_cell inside my_cell);
		ask choosen_cell {
			self.my_cell.is_damaged <- false;
			do die;
		}
		self.differenciation_checked <- true;
		do differentiate("chem1");
	}
	
	action move_up_gradient (string chemical) {
		state <- "follower";
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
	
	
	reflex dispatcher when:(state != "emitter") {
		if (check_chemicals() = false and check_tumor() = false) {
			do random_walk();
		} else if (check_tumor() = true) {
			do kill_damaged_cell();
		} else if (my_cell.chem1 > chem1_threshold) {
			do random_walk();
		} else if (my_cell.chem1 > chem1_sensing_th) {
			if !differenciation_checked {
				bool differentiated <- flip(prob_differentiation);
				if differentiated = true {
					do differentiate("chem2");
				} else {
					do move_up_gradient("chem1");
				}
				self.differenciation_checked <- true;
			} else {
				do move_up_gradient("chem1");
			}
		} else if (my_cell.chem2 > chem2_sensing_th) {
			if !differenciation_checked {
				bool differentiated <- flip(prob_differentiation);
				if differentiated = true {
					do differentiate("chem3");
				} else {
					do move_up_gradient("chem2");
				}
				self.differenciation_checked <- true;
			} else {
				do move_up_gradient("chem2");
			}
		} else if (my_cell.chem3 > chem3_sensing_th) {
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
	float energy <- (rnd(1000) / 1000) * energy_reproduce  update: energy + nutrient_uptake max: energy_reproduce;
	list<tissue_cell> recheable_cells;
	file my_icon <- file("../images/arnold_tuma.png");
	
	init {
		recheable_cells <- my_cell.neighbors where (empty(damaged_cell inside each) and empty(robot inside each));
	}
	
	reflex divide when:!empty(recheable_cells) and (energy >= energy_reproduce) and flip(div_prob) {
		tissue_cell my_cell_tmp <- one_of(recheable_cells);
		ask my_cell_tmp {
			float energy_after_division <- (rnd(1000) / 1000) * myself.energy;
			create damaged_cell number:1 with:[location::self.location, my_cell::tissue_cell(self.location), energy::energy_after_division];
			myself.energy <- myself.energy - energy_after_division;
		}
		
	}
	
	aspect base {
		draw circle(size) color:color;
	}
	
	aspect icon {
		draw my_icon size: 2*size ;
	}
}

experiment tissue_detector type: gui {
	// Inicialization method
	parameter "Inicialization by injection" var:injected init:true category:"Inicialization method";
	
	// Nanorobots parameters
	parameter "Initial number of nanorobots: " var:nb_robots_init init:150 min:1 category:"Nanorobots";
	parameter "Emision rate of chemical 1" var:generation_rt_chem1 init:0.1 min: 0.0 category:"Nanorobots";
	parameter "Emision rate of chemical 2" var:generation_rt_chem2 init:0.1 min: 0.0 category:"Nanorobots";
	parameter "Emision rate of chemical 3" var:generation_rt_chem3 init:0.1 min: 0.0 category:"Nanorobots";
	parameter "Maximum differenciation time" var:max_t_differenciation init:2000.0 category:"Nanorobots";
	
	// Chemical 1 parameters
	parameter "Diffusion rate chemical 1" var:diff_rate_chem1 init:0.3 min: 0.0 max: 1.0 category:"Chemical 1";
	parameter "Evaporation per cycle chemical 1" var:evaporation_per_cycle_chem1 init:0.01 min:0.0 category:"Chemical 1";
	parameter "Chemical 1 threshold" var:chem1_threshold init:0.5 min:0.0 category:"Chemical 1";
	parameter "Sensing chemical 1 threshold" var:chem1_sensing_th init:0.001 category:"Chemical 1";
	
	// Chemical 2 parameters
	parameter "Diffusion rate chemical 2" var:diff_rate_chem2 init:0.3 min: 0.0 max: 1.0 category:"Chemical 2";
	parameter "Evaporation per cycle chemical 2" var:evaporation_per_cycle_chem2 init:0.01 min:0.0 category:"Chemical 2";
	parameter "Sensing chemical 2 threshold" var:chem2_sensing_th init:0.001 category:"Chemical 2";
	
	// Chemical 3 parameters
	parameter "Diffusion rate chemical 3" var:diff_rate_chem3 init:0.3 min: 0.0 max: 1.0 category:"Chemical 3";
	parameter "Evaporation per cycle chemical 3" var:evaporation_per_cycle_chem3 init:0.01 min:0.0 category:"Chemical 3";
	parameter "Sensing chemical 3 threshold" var:chem3_sensing_th init:0.001 category:"Chemical 3";
	
	// Damaged cells parameters
	parameter "Nutrient uptake" var:nutrient_uptake init: 0.1 min:0.01 category:"Damaged cell";
	parameter "Threshold energy to reproduce" var:energy_reproduce init:5.0 category:"Damaged cell";
	
	output {
		display main_display {
			//event mouse_enter action:create_robots;
			grid tissue_cell lines:rgb("grey");
			species robot aspect:base;
			species damaged_cell aspect:base;
		}
		
		display Damaged_cell_information refresh:every(20) {
			chart "Damaged cells evolution" type: series size: {1,1} position: {0, 0} {
				data "number_of_preys" value: nb_damaged_cells color: #red ;
			}
		}
		
		monitor "Number of damaged cells" value: nb_damaged_cells;
	}
	 
}