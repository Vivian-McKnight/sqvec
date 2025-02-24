//! This crate implements the data dynamic array data structure shown in this paper
//! <https://jmmackenzie.io/pdf/mm22-adcs.pdf>
#![allow(dead_code)]
#![allow(internal_features)]
// #![feature(core_intrinsics)]
#![feature(likely_unlikely)]

pub mod sqvec;
pub use sqvec::SqVec;
