/**
* Name: tissuedetectorv2
* Author: san
* Description: 
* Tags: Tag1, Tag2, TagN
*/

model tissuedetectorv2

global torus:false {
	// Initial number of robots
	int nb_robots_init;
	int grid_size <- 50;
	//point center const: true <- { (grid_size / 2),  (grid_size / 2)};
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
	float prob_differentiation <- 0.01;
	
	// Sensing threshold
	float chem1_sensing_th <- 0.00000000000001;
	float chem2_sensing_th <- 0.00000000000001;
	float chem3_sensing_th <- 0.00000000000001;
	
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
			if i < 17 {
				my_cell_ini <- tissue_cell[i+5,j+5];
				create species:robot with:(location:my_cell_ini.location);
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
	
	/*user_command "Create agents here" {
      		create robot number: nb_robots_init {
      			my_cell <- tissue_cell grid_at location::user_location;
      			set location <- my_cell.location;
      		}
   	} */
	
	reflex diffuse {
      diffuse var:chem1 on:tissue_cell proportion: diff_rate_chem1 radius:10 propagation: gradient;
      diffuse var:chem2 on:tissue_cell proportion: diff_rate_chem2 radius:5 propagation: gradient;
      diffuse var:chem3 on:tissue_cell proportion: diff_rate_chem3 radius:5 propagation: gradient;
    }
}

// Agent that represents all the explored tissue (a sample).
// The chemicals will be segregated into its cells.
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
	bool stopped <- false; // A robot will be stopped if it is a guidepost
	bool emitting_chem1 <- false; // True if the robot is emitting the chemical1
	bool emitting_chem2 <- false; // True if the robot is emitting the chemical2
	bool emitting_chem3 <- false; // True if the robot is emitting the chemical3
	
	//tissue_cell my_cell <- one_of (tissue_cell);
	//tissue_cell my_cell <- my_cell_ini;
	tissue_cell my_cell;
	
	init {
		my_cell <- my_cell_ini;
		location <- my_cell.location;
	}
	
	/* AUXILIARY FUNCTIONS */

	bool check_chem1 {
		if my_cell.chem1 > chem1_threshold {
			return true;
		} else {
			return false;
		}
	}
	
	// Check if there is some chemical in the current tissue_cell
	bool check_chemicals {
		tissue_cell my_cell_tmp <-  one_of (my_cell.neighbors where (each.chem1 > chem1_sensing_th and each.chem2 > chem2_sensing_th and each.chem3 > chem1_sensing_th));
		if my_cell_tmp != nil {
			return true;
		} else {
			return false;
		}
	}
	
	// Update chemicals concentrations intop grid cells
	reflex emit_chem1 when:(emitting_chem1=true and stopped){
		tissue_cell(my_cell).chem1 <- tissue_cell(my_cell).chem1+generation_rt_chem1;
	}
	
	reflex emit_chem2 when:(emitting_chem2=true and stopped){
		tissue_cell(my_cell).chem2 <- tissue_cell(my_cell).chem2+generation_rt_chem2;
	}
	
	reflex emit_chem3 when:(emitting_chem3=true and stopped) {
		tissue_cell(my_cell).chem3 <- tissue_cell(my_cell).chem3+generation_rt_chem3;
	}
	
	// Rule 1: If no tumor and no chemical marker -> random walk
	reflex rule1 when:!stopped {
		list<tissue_cell> cells_without_tumor <- my_cell.neighbors where (empty(damaged_cell inside each));
		 if (check_chemicals()=false and (cells_without_tumor != nil)) {
		 	tissue_cell my_cell_tmp <- one_of (shuffle(my_cell.neighbors) where (empty(robot inside each)));
		 	if (my_cell_tmp != nil) {
		 		my_cell <- my_cell_tmp;
		 		location <- my_cell.location;	
		 	} else {
				my_cell <- my_cell;
				location <- my_cell.location;
			} 
		 }
	}
	
	// Rule 2: If detect tumor -> emit chemicar marker 1
	reflex rule2 when:!stopped {
		list<tissue_cell> cells_with_tumor <- my_cell.neighbors where (!(empty(damaged_cell inside each)));
		if (cells_with_tumor != nil) {
			tissue_cell my_cell_tmp <- one_of (cells_with_tumor where (empty(robot inside each)));
			if (my_cell_tmp != nil) {
				damaged_cell tumor_cell <- one_of (damaged_cell inside my_cell_tmp);
				ask tumor_cell {
						my_cell_tmp.is_damaged <- false;
						do die;
					}
					stopped <- true;
					emitting_chem1 <- true;
			} else {
				my_cell <- self.my_cell;
				location <- my_cell.location;
			}
		} else {
			my_cell <- self.my_cell;
			location <- my_cell.location;
		}
	}
	
	
	// Rule 3: If the magnitude of chem1 is greater that threshold -> random walk
	reflex rule3 when:!stopped {
		if (my_cell.chem1 > chem1_threshold) {
			list<tissue_cell> my_cells_tmp <- (my_cell.neighbors where (empty(robot inside each)));
			if (my_cells_tmp != nil) {
				my_cell <- one_of (my_cells_tmp);
				location <- my_cell.location;
			} else {
				my_cell <- self.my_cell;
				location <- my_cell.location;
			}
			//my_cell <- one_of (my_cell.neighbors);
		 	//self.location <- my_cell.location;
		}
	} 
	
	// Rule 4: if chem1 is detected -> move up the gradient and differentiate with prob P.
	reflex rule4 when:!stopped {
		if (my_cell.chem1 > 0.0) {
			bool differentiated <- flip(prob_differentiation);
			if differentiated=true {
				stopped <- true;
				emitting_chem2 <- true;
			} else {
				list<tissue_cell> my_cells_tmp <- (my_cell.neighbors) where (empty(robot inside each));
				if my_cells_tmp != nil {
					tissue_cell my_cell_maximize <- (my_cells_tmp) with_max_of (each.chem1);
					if (my_cell_maximize != nil) {
						my_cell <- my_cell_maximize;
						location <- my_cell.location;
					} else {
						my_cell <- self.my_cell;
						location <- my_cell.location;
					}
				} else {
					my_cell <- self.my_cell;
					location <- my_cell.location;
				}
			}
		}
	} 
	
	reflex rule5 when:!stopped {
		if (my_cell.chem2 > 0.0) {
			bool differentiated <- flip(prob_differentiation);
			if differentiated=true {
				stopped <- true;
				emitting_chem3 <- true;
			} else {
				list<tissue_cell> my_cells_tmp <- (my_cell.neighbors) where (empty(robot inside each));
				if my_cells_tmp != nil {
					tissue_cell my_cell_maximize <- (my_cells_tmp) with_max_of (each.chem2);
					if (my_cell_maximize != nil) {
						my_cell <- my_cell_maximize;
						location <- my_cell.location;
					} else {
						my_cell <- self.my_cell;
						location <- my_cell.location;
					}
				} else {
					my_cell <- self.my_cell;
					location <- my_cell.location;
				}
			}
		}
	}
	
	reflex rule6 when:!stopped {
		if (my_cell.chem3 > 0.0) {
			list<tissue_cell> my_cells_tmp <- (my_cell.neighbors) where (empty(robot inside each));
			if my_cells_tmp != nil {
				tissue_cell my_cell_maximize <- (my_cells_tmp) with_max_of (each.chem2);
				if (my_cell_maximize != nil) {
					self.my_cell <- my_cell_maximize;
					self.location <- my_cell.location;
				} else {
					my_cell <- self.my_cell;
					location <- my_cell.location;
				}
			} else {
				my_cell <- self.my_cell;
				location <- my_cell.location;
			}
			
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
	
	reflex divide when:(recheable_cells != nil) and (flip(div_prob)) {
		tissue_cell my_cell_tmp <- one_of(recheable_cells);
		ask my_cell_tmp {
			create damaged_cell number:1 with:[location::self.location, my_cell::tissue_cell(self.location)];
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
	parameter "Initial number of nanorobots: " var:nb_robots_init init:270 min:1 category:"Nanorobots";
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