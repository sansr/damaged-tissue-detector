/**
 *  new
 *  Author: Sandra Saez
 *  Description: 
 */

model tissue_detector

global {
	// Initial number of robots
	int nb_robots_init;
	int grid_size <- 50;
	
	// Diffusion  rate from input
	float diff_rate_chem1;
	float diff_rate_chem2;
	float diff_rate_chem3;
	
	// Emision rate of chemicals.
	float generation_rt_chem1;
	float generation_rt_chem2;
	float generation_rt_chem3;
	
	// Read image of the tissue
	file map_init <- image_file("../images/damaged_tissue2.png");
	
	// Set up the world
	init {
		matrix init_data <- map_init as_matrix {grid_size,grid_size};
		create robot number:nb_robots_init;
		ask tissue_cell {
			loop cell over:init_data {
				rgb cell_color <- rgb(init_data[grid_x,grid_y]);
				if cell_color != #white {
					self.is_damaged <- true;
				}
			}
			
			// Create the damaged_cell agents
			if self.is_damaged {
				create damaged_cell number:1 {
					self.location <- myself.location;
				}
			}
		}	
	}
	
	reflex diffuse {
      diffuse var:chem1 on:tissue_cell proportion: 0.9 radius:2 propagation: gradient;
   }
}

// Agent that represents all the explored tissue (a sample).
// The chemicals will be segregated into its cells.
grid tissue_cell height:grid_size width:grid_size neighbors:8 {
	bool is_damaged <- false;
	rgb color <- rgb(int(255 * (1 - chem1)), 255 * (1 - chem2), int(255 * (1 - chem3))) update:rgb(int(255 * (1 - chem1)), 255 * (1 - chem2), int(255 *(1 - chem3))) ;
	
	// Concentrations of the three different chemicals.
	float chem1 <- 0.0; 
	float chem2 <- 0.0;
	float chem3 <- 0.0;
	
	// List of neighbours at distance 1
	list<tissue_cell> neighbors <- self neighbors_at 1;
}

// Agent that represents each of the nanorobots.
species robot {
	float size <- 0.5;
	rgb color <- #blue;
	bool stopped <- false; // A robot will be stopped if it is a guidepost
	bool emitting_chem1 <- false; // True if the robot is emitting the chemical1
	bool emitting_chem2 <- false; // True if the robot is emitting the chemical2
	bool emitting_chem3 <- false; // True if the robot is emitting the chemical3
	
	tissue_cell my_cell <- one_of (tissue_cell);
	
	init {
		location <- my_cell.location;
	}
	
	/* AUXILIARY FUNCTIONS */
	
	// Choose a cell without robots. Avoid overlaps.
	tissue_cell choose_cell_without_robot {
		tissue_cell my_cell_tmp <- shuffle(my_cell.neighbors) first_with (!(empty (robot inside (each))));
		if my_cell_tmp = nil {
			return one_of (my_cell.neighbors);
		} else {
			list<tissue_cell> cells <- my_cell.neighbors where (empty(robot inside (each)));
			if cells != nil {
				return one_of (cells);
			} else {
				return my_cell;
			}
		}
	}
	
	// Choose one cell to move
	tissue_cell choose_cell {
		tissue_cell my_cell_tmp <- shuffle(my_cell.neighbors) first_with (!(empty (damaged_cell inside (each))));
		if my_cell_tmp != nil {
			if empty(robot inside my_cell_tmp) = true {
				return my_cell_tmp;
			} else {
				tissue_cell cell_without_robot <- choose_cell_without_robot();
				return cell_without_robot;
			}	
		} else {
			tissue_cell cell_without_robot <- choose_cell_without_robot();
			return cell_without_robot;
		}
	}
	
	// Check if there is some chemical in the current tissue_cell
	action check_chemicals type:bool {
		if my_cell.chem1 != 0.0 or my_cell.chem2 != 0.0 or my_cell.chem3 != 0.0{
			return true;
		} else {
			return false;
		}
	}
	
	reflex emit_chem1 when:emitting_chem1=true{
		tissue_cell(location).chem1 <- tissue_cell(location).chem1+generation_rt_chem1;
	}
	
	reflex emit_chem2 when:emitting_chem2=true{
		tissue_cell(location).chem2 <- tissue_cell(location).chem2+generation_rt_chem2;
	}
	
	reflex emit_chem3 when:emitting_chem3=true{
		tissue_cell(location).chem3 <- tissue_cell(location).chem3+generation_rt_chem3;
	}
	
	/* EAT TUMOR RULE SET */
	
	// Rule 1: If no tumor and no chemical marker -> random walk
	reflex move when:!stopped {
		my_cell <- choose_cell();
		location <- my_cell.location; 
	}
	
	// Rule 2: If detect tumor -> emit chemicar marker 1
	reflex detect_tumor when:!stopped {
		// If a there exists a damaged cell inside the robot cell -> kill it
		list<damaged_cell> cell_to_kill <- damaged_cell inside my_cell;
		if empty(cell_to_kill) = false {
			stopped <- true;
			ask one_of (cell_to_kill) {
				do die;
			}
			emitting_chem1 <- true;
		}
	}
	
	aspect base {
		draw circle(size) color:color;
	}
}

species damaged_cell {
	aspect base {
		draw circle(0.5) color:#red;
	}
}

experiment tissue_detector type: gui {
	parameter "Initial number of nanorobots: " var:nb_robots_init init:20 min:1 category:"Nanorobots";
	parameter "Emision rate of chemical 1" var:generation_rt_chem1 init:0.1 min: 0.0 category:"Nanorobots";
	parameter "Emision rate of chemical 2" var:generation_rt_chem2 init:0.1 min: 0.0 category:"Nanorobots";
	parameter "Emision rate of chemical 3" var:generation_rt_chem3 init:0.1 min: 0.0 category:"Nanorobots";
	parameter "Diffusion rate of chemical 1" var:diff_rate_chem1 init:0.5 min: 0.0 max: 1.0 category:"Nanorobots";
	parameter "Diffusion rate of chemical 2" var:diff_rate_chem2 init:0.5 min: 0.0 max: 1.0 category:"Nanorobots";
	parameter "Diffusion rate of chemical 3" var:diff_rate_chem3 init:0.5 min: 0.0 max: 1.0 category:"Nanorobots";
	output {
		display main_display {
			grid tissue_cell lines:rgb("grey");
			species robot aspect:base;
			species damaged_cell aspect:base;
		}
	}
}
