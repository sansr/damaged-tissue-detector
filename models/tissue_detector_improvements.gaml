/**
* Name: tissue_detector_improvements
* Author: Sandra Saez
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
	float max_lost_timer;
	
	// Grid global variables
	int grid_size;
	
	// Damaged cell global variables
	float energy_reproduce;
	float nutrient_uptake;
	int nb_damaged_cells -> {length (damaged_cell)};
	int nb_robots -> {length (robot)};
	int total_differentiated;
	
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
	float div_prob;
	
	// Read image of the tissue
	string relative_path <- "../images/damaged_tissue2.png";
	file map_init <- image_file(relative_path);
	
	// Set up the world
	init {
		matrix init_data <- map_init as_matrix {grid_size,grid_size};
		
		if injected = false {
			int j <- 0; // row
			int i <- 0; // column
			int side_dimension <- int(sqrt(world.nb_robots_init));
			loop times:nb_robots_init+1 {
				if i < side_dimension {
					my_cell_ini <- tissue_cell[i,j];
					create robot with:[my_cell::my_cell_ini, location::my_cell_ini.location, my_first_cell::my_cell_ini];
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
	rgb color <- rgb(int(255 * (1 - chem1)), 255 * (1 - chem2), int(255 * (1 - chem3))) update:rgb(int(255 * (1 - chem1)), 255 * (1 - chem2), int(255 *(1 - chem3)));
	
	// Concentrations of the three different chemicals.
	float chem1 <- 0.0 update: (chem1<=evaporation_per_cycle_chem1) ? 0.0 : chem1-evaporation_per_cycle_chem1; 
	float chem2 <- 0.0 update: (chem2<=evaporation_per_cycle_chem2) ? 0.0 : chem2-evaporation_per_cycle_chem2;
	float chem3 <- 0.0 update: (chem3<=evaporation_per_cycle_chem3) ? 0.0 : chem3-evaporation_per_cycle_chem3;
	
	// List of neighbours at distance 1
	list<tissue_cell> neighbors <- self neighbors_at 1;
	
	action create_robots {
		create robot number:nb_robots_init with: [location::self.location, my_cell::tissue_cell(location), my_first_cell::tissue_cell(location)];
	}
	
	action create_robots_from_cell {
		int i <- 0;
		int j <- 0;
		int side_dimension <- int(sqrt(world.nb_robots_init));
		int iter <- 0;
		loop times:side_dimension {
			tissue_cell cell_tmp;
			loop times:side_dimension {
				if i < side_dimension {
					cell_tmp <- tissue_cell[grid_x+i, grid_y+j];
					create robot with:[my_cell::cell_tmp, location::cell_tmp.location, my_first_cell::cell_tmp];
					i <- i+1;
				}
			}
			i <- 0;
			j <- j + 1;
		}
	} 
	
	user_command create_robots_here action:create_robots;
	user_command create_robots_from_here action:create_robots_from_cell;
}

// Agent that represents each of the nanorobots.
species robot {
	float size <- 0.5;
	rgb color <- #lime;
	string chemicals <- "chem1" among: ["chem1", "chem2", "chem3"];
	string emitting <- "NO" among: ["NO", "CHEM1", "CHEM2", "CHEM3"];
	string state <- "wander" among: ["wander", "follower", "emitter", "going_injection_site"];
	bool differenciation_checked <- false;
	
	// Temporal differenciation
	float max_time <- max_t_differenciation;
	float local_timer <- 0.0 max: max_time;
	
	// Loose timer
	float lost_timer <- max_lost_timer;
	
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
	
	/*************************************** */
	/*           LIST OF BEHAVIOURS            */
	/*************************************** */
	
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
		total_differentiated <- total_differentiated-1;
		self.local_timer <- 0.0; // Reset the timer
		self.differenciation_checked <- false; // To have the opportunity to differenciate again
	}
	
	reflex i_am_lost when:lost_timer=0 {
		do die;
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
		total_differentiated <- total_differentiated+1; 
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
	/*          EAT TUMOR RULE SET           */
	/*************************************** */
	
	
	reflex eat_tumor_rule_set when:(state != "emitter" and state != "going_injection_site") {
		if (check_chemicals() = false and check_tumor() = false) {
			lost_timer <- lost_timer - 1.0;
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
	
	/*************************************** */
	/*     GOING INJECTION SITE RULE     */
	/*************************************** */
	
	reflex go_back_injection_site_rule_set when: (nb_damaged_cells = 0) {
		emitting <- "NO";
		state <- "going_injection_site";

		int xoff <- my_first_cell.grid_x - my_cell.grid_x;
		int yoff <- my_first_cell.grid_y - my_cell.grid_y;
		
		int newx;
		if xoff = 0 {
			newx <- my_cell.grid_x;
		} else {
			newx <- my_cell.grid_x + round(xoff/abs(xoff));
		}
	
		int newy;
		if yoff = 0 {
			newy <- my_cell.grid_y;
		} else {
			newy <- my_cell.grid_y + round(yoff/abs(yoff));
		}
		
		my_cell <- tissue_cell[newx, newy];
		location <- my_cell.location;
		//do goto target:my_first_cell;
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
}

experiment tissue_detector type: gui {
	//Grid
	parameter "Grid dimensions" var:grid_size init:50 category:"Grid";
	parameter "Image" var:relative_path init:"../images/damaged_tissue2.png" category:"Grid";
	
	// Inicialization method
	parameter "Inicialization by injection" var:injected init:true category:"Inicialization method";
	
	// Nanorobots parameters
	parameter "Initial number of nanorobots: " var:nb_robots_init init:289 min:1 category:"Nanorobots";
	parameter "Differenciation probability" var:prob_differentiation init:0.01 category:"Nanorobots";
	parameter "Emision rate of chemical 1" var:generation_rt_chem1 init:0.1 min: 0.0 category:"Nanorobots";
	parameter "Emision rate of chemical 2" var:generation_rt_chem2 init:0.1 min: 0.0 category:"Nanorobots";
	parameter "Emision rate of chemical 3" var:generation_rt_chem3 init:0.1 min: 0.0 category:"Nanorobots";
	parameter "Maximum differenciation time" var:max_t_differenciation init:500.0 category:"Nanorobots";
	parameter "Lost timer" var:max_lost_timer init:4000.0 category:"Nanorobots";
	
	// Chemical 1 parameters
	parameter "Diffusion rate chemical 1" var:diff_rate_chem1 init:0.3 min: 0.0 max: 1.0 category:"Chemical 1";
	parameter "Degradation per cycle chemical 1" var:evaporation_per_cycle_chem1 init:0.01 min:0.0 category:"Chemical 1";
	parameter "Chemical 1 threshold" var:chem1_threshold init:0.5 min:0.0 category:"Chemical 1";
	parameter "Sensing chemical 1 threshold" var:chem1_sensing_th init:0.001 category:"Chemical 1";
	
	// Chemical 2 parameters
	parameter "Diffusion rate chemical 2" var:diff_rate_chem2 init:0.3 min: 0.0 max: 1.0 category:"Chemical 2";
	parameter "Degradation per cycle chemical 2" var:evaporation_per_cycle_chem2 init:0.01 min:0.0 category:"Chemical 2";
	parameter "Sensing chemical 2 threshold" var:chem2_sensing_th init:0.001 category:"Chemical 2";
	
	// Chemical 3 parameters
	parameter "Diffusion rate chemical 3" var:diff_rate_chem3 init:0.3 min: 0.0 max: 1.0 category:"Chemical 3";
	parameter "Degradation per cycle chemical 3" var:evaporation_per_cycle_chem3 init:0.01 min:0.0 category:"Chemical 3";
	parameter "Sensing chemical 3 threshold" var:chem3_sensing_th init:0.001 category:"Chemical 3";
	
	// Damaged cells parameters
	parameter "Nutrient uptake" var:nutrient_uptake init: 0.1 min:0.01 category:"Damaged cell";
	parameter "Threshold energy to reproduce" var:energy_reproduce init:5.0 category:"Damaged cell";
	parameter "Division probability" var:div_prob init:0.001 category:"Damaged cell";
	
	output {
		display main_display {
			//event mouse_enter action:create_robots;
			grid tissue_cell lines:rgb("grey");
			species robot aspect:base;
			species damaged_cell aspect:base;
		}
		
		display Damaged_cell_information refresh:every(20) {
			chart "Damaged cells evolution" type: series size: {1,1} position: {0, 0} {
				data "number_of_damaged_cells" value: nb_damaged_cells color: #red ;
			}
		}
		display Robots_differentiated refresh:every(20) {
			chart "Robots differenciated" type: series size: {1,1} position: {0, 0} {
				data "robots_differentiated" value:total_differentiated color: #purple ;
			}
		}
		display Robots_information refresh:every(20) {
			chart "Robots evolution" type: series size: {1,1} position: {0, 0} {
				data "number_of_robots" value: nb_robots color: #blue ;
			}
		}
		
		monitor "Number of damaged cells" value: nb_damaged_cells;
		monitor "Number of robots" value: nb_robots;
	}
	 
}