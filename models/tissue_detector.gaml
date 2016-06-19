/**
 *  new
 *  Author: Sandra Saez
 *  Description: 
 */

model tissue_detector

global {
	// Initial number of robots
	int nb_robots_init;
	
	// Emission rate from input
	float rt_chem1 <- 0.0;
	float rt_chem2 <- 0.0;
	float rt_chem3 <- 0.0;
	
	// Read image of the tissue
	file map_init <- image_file("../images/damaged_tissue2.png");
	
	// Set up the world
	init {
		matrix init_data <- map_init as_matrix {50,50};
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
					set location <- myself.location;
				}
			}
		}	
	}
}

// Agent that represents all the explored tissue (a sample).
// The chemicals will be segregated into its cells.
grid tissue_cell height:50 width:50 neighbours:8 {
	bool is_damaged <- false;
	rgb color <- rgb(int(255 * (1 - chem1)), 255 * (1 - chem2), int(255 * (1 - chem3))) update:rgb(int(255 * (1 - chem1)), 255 * (1 - chem2), int(255 *(1 - chem3))) ;
	
	// Emision rate of chemicals.
	float generation_rt_chem1 <- 0.0;
	float generation_rt_chem2 <- 0.0;
	float generation_rt_chem3 <- 0.0;
	
	// Concentrations of the three different chemicals.
	//list<robot> local_robot <- one_of (robot inside {grid_x, grid_y});
	float chem1 <- 0.0 update:chem1+generation_rt_chem1;
	float chem2 <- 0.0 update:chem2+generation_rt_chem2;
	float chem3 <- 0.0 update:chem3+generation_rt_chem3;
	
	// List of neighbours at distance 1
	list<tissue_cell> neighbours <- self neighbours_at 1;
}

// Agent that represents each of the nanorobots.
species robot skills:[moving]{
	float size <- 0.5;
	rgb color <- #blue;
	bool stopped <- false; // A robot will be stopped if it is a guidepost
	bool emitting <- false; // True if the robot is emitting some chemical
	
	tissue_cell my_cell <- one_of (tissue_cell);
	list<damaged_cell> reachable_damaged_cells update:damaged_cell inside (my_cell.neighbours);
	
	init {
		location <- my_cell.location;
	}
	
	/* AUXILIARY FUNCTIONS */
	
	// Check if there is some chemical in the current tissue_cell
	action check_chemicals type:bool {
		if my_cell.chem1 != 0.0 or my_cell.chem2 != 0.0 or my_cell.chem3 != 0.0{
			return true;
		} else {
			return false;
		}
	}
	
	// Chemicals diffusion
	/*reflex diff when:emitting {
		list cells_where_diffuse <- tissue_cell where (each.grid_x < 2);
		loop n over:cells_where_diffuse {
			n.color <- #yellow;
		}      
    }*/
	
	/* EAT TUMOR RULE SET */
	
	// Rule 1: If no tumor and no chemical marker -> random walk
	reflex move when:!stopped {
		if empty(reachable_damaged_cells) and !check_chemicals() {
			my_cell <- one_of (my_cell.neighbours);
			location <- my_cell.location;
		}
	}
	
	// Rule 2: If detect tumor -> emit chemicar marker 1
	reflex detect_tumor when:!stopped {
		// If a neighbour cell is damaged
		if !empty(reachable_damaged_cells) {
			ask one_of (reachable_damaged_cells) {
				do die;
			}
					
			// Stop and start emiting chem1
			stopped <- true;
			emitting <- true;
			my_cell.generation_rt_chem1 <- 0.01;
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
	parameter "Generation ratio of chemical 1" var:rt_chem1 init:0.01 category:"Nanorobots";
	parameter "Generation ratio of chemical 2" var:rt_chem2 init:0.01 category:"Nanorobots";
	parameter "Generation ratio of chemical 3" var:rt_chem3 init:0.01 category:"Nanorobots";
	output {
		display main_display {
			grid tissue_cell lines:rgb("grey");
			species robot aspect:base;
			species damaged_cell aspect:base;
		}
	}
}
