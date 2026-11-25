//! Command-line builder for coarse-grained iSoLF membranes, vesicles, tubules,
//! AICG2+ proteins, and 1D bilayer umbrella-sampling campaigns for GENESIS.
//!
//! ```sh
//! isolf membrane --lipids "POPC=3,POPS=1" --count 512 --inp --vmd -o ./sim
//! isolf membrane --lipids POPC --box 20 -p channel.pdb --lipid-protein lp.txt --inp -o ./prot_mem
//! isolf vesicle  --lipids POPC --radius 20 --top -o ./ves
//! isolf tubule   --lipids "DPPC=1,POPC=3" --phase raft --radius 10 --length 25 -o ./tube
//! isolf protein  1aki.pdb --inp --vmd -o ./prot
//! isolf umbrella setup --probe TRP -L POPC --z 0.0..4.0:0.1 -o ./us_trp
//! isolf umbrella cv --discard 50000 -o ./us_trp
//! isolf params   --params trial.txt -o ./us_trp
//! isolf lipids
//! ```

mod anim;
mod cli;
mod report;
mod runlog;
mod update;

fn main() {
    cli::main();
}
