use std::{io::{Read, Write}, ops::Range};
use std::sync::atomic::{AtomicU16, AtomicU32, AtomicU64, AtomicU8};
use std::sync::atomic::Ordering::Relaxed;

#[derive(PartialEq, Eq, Debug, Copy, Clone)]
pub enum ColorVecValue {
    Single(usize),
    Multiple,
    None,
} 

pub trait ColorStorage {
    fn get_color(&self, colex: usize) -> ColorVecValue;
    fn get_color_of_range(&self, range: Range<usize>) -> ColorVecValue;
}

pub trait MySerialize {
    fn serialize(&self, out: &mut impl Write);
    fn load(input: &mut impl Read) -> Box<Self>;
}

pub trait AtomicColorVec{
    // Represents None as the max value of the atomic type

    fn update(&self, i: usize, x: usize); 
    fn read(&self, i: usize) -> ColorVecValue;
    fn new(len: usize) -> Self; // Stores a None (=max_value()) to each position
}

pub trait AtomicUint {
    fn load(&self, order: std::sync::atomic::Ordering) -> usize;
    fn fetch_update<F: Fn(usize) -> usize>(&self, f: F);
    fn max_value() -> usize;
    fn new(val: usize) -> Self;
}

impl AtomicUint for AtomicU8 {
    fn load(&self, order: std::sync::atomic::Ordering) -> usize {
        self.load(order) as usize
    }
    fn fetch_update<F: Fn(usize) -> usize>(&self, f: F) {
        self.fetch_update(Relaxed, Relaxed, |x| Some(f(x as usize) as u8)).unwrap();
    }
    fn max_value() -> usize {
        u8::MAX as usize
    }
    fn new(val: usize) -> Self {
        AtomicU8::new(val as u8)
    }
}

impl AtomicUint for AtomicU16 {
    fn load(&self, order: std::sync::atomic::Ordering) -> usize {
        self.load(order) as usize
    }
    fn fetch_update<F: Fn(usize) -> usize>(&self, f: F) {
        self.fetch_update(Relaxed, Relaxed, |x| Some(f(x as usize) as u16)).unwrap();
    }
    fn max_value() -> usize {
        u16::MAX as usize
    }
    fn new(val: usize) -> Self {
        AtomicU16::new(val as u16)
    }
}

impl AtomicUint for AtomicU32 {
    fn load(&self, order: std::sync::atomic::Ordering) -> usize {
        self.load(order) as usize
    }
    fn fetch_update<F: Fn(usize) -> usize>(&self, f: F) {
        self.fetch_update(Relaxed, Relaxed, |x| Some(f(x as usize) as u32)).unwrap();
    }
    fn max_value() -> usize {
        u32::MAX as usize
    }
    fn new(val: usize) -> Self {
        AtomicU32::new(val as u32)
    }
}

impl AtomicUint for AtomicU64 {
    fn load(&self, order: std::sync::atomic::Ordering) -> usize {
        self.load(order) as usize
    }
    fn fetch_update<F: Fn(usize) -> usize>(&self, f: F) {
        self.fetch_update(Relaxed, Relaxed, |x| Some(f(x as usize) as u64)).unwrap();
    }
    fn max_value() -> usize {
        u64::MAX as usize
    }
    fn new(val: usize) -> Self {
        AtomicU64::new(val as u64)
    }
}