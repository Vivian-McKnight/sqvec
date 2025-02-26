use std::alloc::{self, Layout, alloc, dealloc, realloc};
use std::fmt::Debug;
use std::marker::PhantomData;
use std::ops::{Index, IndexMut};
use std::ptr::NonNull;
use std::{hint, mem};

// todo
// - add shrinking the array to size
// - efficient item drop
// - c interop
// - size checks/debug assertions
// - initialisation macro
// - clone trait
// - from slice/box trait
// - intoiter
// - insert
// - remove
// - store seglen instead of log_seglen?
// - zero sized types
// - benchmark different iteration strategies
// - benchmark random operations against vec
// - usize as index?

/// A non-contiguos extensible array that achieves asymptotically optimal
/// space overhead of O(√n) while maintaining comparable time efficiency**.
/// Specification: <https://jmmackenzie.io/pdf/mm22-adcs.pdf>
pub struct SqVec<T> {
    dope: NonNull<NonNull<T>>,
    dope_cap: u32,
    len: u32,
    alloc_seg_count: u32,
    marker: PhantomData<T>,
}

#[macro_export]
macro_rules! sqvec {
    ( $( $x:expr ),* ) => {
        {
            let mut temp_vec = SqVec::new();
            $(
                temp_vec.push($x);
            )*
            temp_vec
        }
    };
}

impl<T> SqVec<T> {
    /// Constructs a new, empty `SqVec<T>`.
    ///
    /// The sqvec will not allocate until elements are pushed onto it.
    pub fn new() -> Self {
        assert!(
            std::mem::size_of::<T>() != 0,
            "This data structure does not support ZST's"
        );
        Self {
            dope: NonNull::dangling(),
            dope_cap: 0,
            len: 0,
            alloc_seg_count: 0,
            marker: PhantomData,
        }
    }

    /// calculates the dope_cap from `self.alloc_seg_count`
    fn dope_cap(&self) -> u32 {
        todo!()
    }

    /// grows the dope vector to the correct size to fit the next set of segments
    fn grow_dope(&mut self) {
        let (new_cap, new_layout) = if self.dope_cap == 0 {
            (2, unsafe {
                Layout::array::<NonNull<T>>(2).unwrap_unchecked()
            })
        } else {
            let new_cap = self.dope_cap + Self::segs_of_len(Self::lsln(self.len - 1));
            let new_layout =
                Layout::array::<NonNull<T>>(new_cap as usize).expect("Allocation too large");
            (new_cap, new_layout)
        };

        let handle = if self.dope_cap == 0 {
            unsafe { alloc(new_layout) }
        } else {
            // this layout cannot error because its size is based on an allocation that already exists
            let old_layout =
                unsafe { Layout::array::<NonNull<T>>(self.dope_cap as usize).unwrap_unchecked() };
            let old_handle = self.dope.as_ptr() as *mut u8;
            unsafe { realloc(old_handle, old_layout, new_layout.size()) }
        };

        self.dope = NonNull::new(handle as *mut NonNull<T>)
            .unwrap_or_else(|| alloc::handle_alloc_error(new_layout));
        self.dope_cap = new_cap;
    }

    /// Allocates a new segment of length `1 << self.log_seglen` at index `self.alloc_seg_count` in dope.
    fn alloc_seg(&mut self) {
        let layout = Layout::array::<T>(1 << Self::lsln(self.len)).expect("Allocation too large");
        let handle = unsafe { alloc(layout) };
        let ptr =
            NonNull::new(handle as *mut T).unwrap_or_else(|| alloc::handle_alloc_error(layout));
        unsafe { self.dope.add(self.alloc_seg_count as usize).write(ptr) };
        self.alloc_seg_count += 1;
    }

    /// Adds `val` to the end of the SqVec, increasing its length by 1.
    pub fn push(&mut self, val: T) {
        let (segnum, offset) = Self::mapping(self.len);
        if offset == 0 {
            if segnum == self.alloc_seg_count {
                if segnum == self.dope_cap {
                    self.grow_dope();
                }
                self.alloc_seg();
            }
        }
        unsafe { self.item_ptr_inner(segnum, offset).write(val) };
        self.len += 1;
    }

    /// Removes the last value on the SqVec and returns it in an Option.
    /// If the SqVec is empty, Returns None
    pub fn pop(&mut self) -> Option<T> {
        if self.len == 0 {
            None
        } else {
            self.len -= 1;
            Some(unsafe { self.item_ptr(self.len).read() })
        }
    }

    /// swaps the elements at index's ix1 and ix2
    pub fn swap(&mut self, ix1: u32, ix2: u32) {
        assert!(ix1 < self.len && ix2 < self.len, "Index out of range");
        unsafe { self.swap_unchecked(ix1, ix2) };
    }

    /// swaps the elements at index's ix1 and ix2 without checking ix1 or ix2 in range
    pub unsafe fn swap_unchecked(&mut self, ix1: u32, ix2: u32) {
        unsafe { mem::swap(self.item_ptr(ix1).as_mut(), self.item_ptr(ix2).as_mut()) };
    }

    /// Returns the length of the SqVec (number of elements)
    #[inline]
    pub fn len(&self) -> usize {
        self.len as usize
    }

    /// Returns true if the sqvec is empty and false otherwise
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Maps the given index to the segment number and index within the
    /// segment (offset) where the element at index `v` can be found.
    #[inline]
    fn mapping(v: u32) -> (u32, u32) {
        let b = if v == 0 {
            1
        } else {
            (33 - v.leading_zeros()) >> 1
        };
        let segnum = (v >> b) + (1 << (b - 1)) - 1;
        let offset = v & ((1 << b) - 1);
        (segnum, offset)
    }

    /// Gives the dope index of the first segment of length (1 << log_seglen)
    #[inline]
    fn first_seg(log_seglen: u32) -> u32 {
        if log_seglen == 1 {
            0
        } else {
            3 * (1 << (log_seglen - 2)) - 1
        }
    }

    /// Gives the dope index of the last segment of length (1 << log_seglen)
    #[inline]
    fn last_seg(log_seglen: u32) -> u32 {
        3 * (1 << (log_seglen - 1)) - 2
    }

    /// Returns the max number of segments of lenght (1 << log_seglen) that can
    /// be allocated
    #[inline]
    fn segs_of_len(log_seglen: u32) -> u32 {
        if log_seglen == 1 {
            2
        } else {
            3 * (1 << (log_seglen - 2))
        }
    }

    /// Returns the pointer to the item in the SqVec at index `ix` (Does not check if ix is < self.len)
    #[inline]
    unsafe fn item_ptr(&self, ix: u32) -> NonNull<T> {
        let (segnum, offset) = Self::mapping(ix);
        unsafe { self.item_ptr_inner(segnum, offset) }
    }

    #[inline]
    unsafe fn item_ptr_inner(&self, segnum: u32, offset: u32) -> NonNull<T> {
        unsafe { self.dope.add(segnum as usize).read().add(offset as usize) }
    }

    pub fn iter(&self) -> Iter<T> {
        Iter { data: &self, ix: 0 }
    }

    pub fn iter_mut(&mut self) -> IterMut<T> {
        IterMut { data: self, ix: 0 }
    }

    /// calculates the log base 2 segment length of the
    /// segment that the element at index `v` would be located
    #[inline]
    fn lsln(v: u32) -> u32 {
        if v == 0 {
            return 1;
        }
        (32 - v.leading_zeros()).div_ceil(2)
    }

    /// clears the sqvec dropping all present elements but
    /// maintaining the currently allocated space
    fn clear(&mut self) {
        if mem::needs_drop::<T>() {
            for i in 0..self.len {
                unsafe { self.item_ptr(i).drop_in_place() };
            }
        }
        self.len = 0;
    }

    pub fn iter_t(&self) -> IterT<T> {
        let (end_seg_ix, end_element_ix) = Self::mapping(self.len - 1);
        IterT {
            dope_iter: unsafe {
                core::slice::from_raw_parts(
                    self.dope.as_ptr() as *const NonNull<T>,
                    end_seg_ix as usize,
                )
            }
            .iter()
            .copied(),
            seg_iter: unsafe {
                core::slice::from_raw_parts(
                    self.dope.read().as_ptr() as *const T,
                    2.min(self.len()),
                )
            }
            .iter(),
            end_seg_ptr: unsafe { self.dope.add(end_seg_ix as usize).read() },
            end_seg_len: end_element_ix + 1,
            segs_left: 2,
            seglen: 2,
            marker: PhantomData,
        }
    }
}

#[derive(Debug)]
pub struct IterT<'a, T: 'a> {
    dope_iter: core::iter::Copied<core::slice::Iter<'a, NonNull<T>>>,
    seg_iter: core::slice::Iter<'a, T>,
    end_seg_ptr: NonNull<T>,
    end_seg_len: u32,
    /// number of segments left of length seglen
    segs_left: u32,
    seglen: u32,
    marker: PhantomData<&'a T>,
}

impl<'a, T> Iterator for IterT<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(item) = self.seg_iter.next() {
            Some(item)
        } else {
            if let Some(ptr) = self.dope_iter.next() {
                if self.segs_left == 0 {
                    self.seglen <<= 1;
                    self.segs_left = (3 * self.seglen) >> 1; // self.seglen will only ever be >= 4 at this point 👍
                }

                if ptr == self.end_seg_ptr {
                    self.seg_iter = unsafe {
                        core::slice::from_raw_parts(
                            ptr.as_ptr() as *const T,
                            self.end_seg_len as usize,
                        )
                    }
                    .iter();
                } else {
                    self.seg_iter = unsafe {
                        core::slice::from_raw_parts(ptr.as_ptr() as *const T, self.seglen as usize)
                    }
                    .iter();
                }

                self.segs_left -= 1;
                self.seg_iter.next()
            } else {
                None
            }
        }
    }
}

struct DopeIter<'a, T: 'a> {
    ptr: NonNull<NonNull<T>>,
    end: NonNull<NonNull<T>>,
    marker: PhantomData<&'a T>,
}

impl<'a, T> DopeIter<'a, T> {
    // fn new()
}

impl<'a, T> Iterator for DopeIter<'a, T> {
    type Item = NonNull<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.ptr > self.end {
            None
        } else {
            let out = unsafe { self.ptr.read() };
            self.ptr = unsafe { self.ptr.add(1) };
            Some(out)
        }
    }
}

impl<T> Drop for SqVec<T> {
    fn drop(&mut self) {
        if self.dope_cap == 0 {
            return;
        }
        // drop items in SqVec
        if mem::needs_drop::<T>() {
            for i in 0..self.len {
                unsafe { self.item_ptr(i).drop_in_place() };
            }
        }

        let dope_slice = unsafe {
            std::slice::from_raw_parts(self.dope.as_ptr(), self.alloc_seg_count as usize)
        };
        let mut ix1 = 0_usize;
        let max_lsl = u32::ilog2((self.dope_cap + 2) / 3) + 1;
        for lsl in 1..=max_lsl {
            let ix2 = Self::last_seg(lsl).min(self.alloc_seg_count.saturating_sub(1)) as usize;
            let layout = Layout::array::<T>(1 << lsl).unwrap();
            for &ptr in dope_slice[ix1..=ix2].iter() {
                unsafe { dealloc(ptr.as_ptr() as *mut u8, layout) };
            }
            ix1 = ix2 + 1;
        }

        unsafe {
            dealloc(
                self.dope.as_ptr() as *mut u8,
                Layout::array::<NonNull<T>>(self.dope_cap as usize).unwrap(),
            )
        };
    }
}

pub struct IterMut<'a, T> {
    data: &'a SqVec<T>,
    ix: u32,
}

impl<'a, T> Iterator for IterMut<'a, T> {
    type Item = &'a mut T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.ix == self.data.len() as u32 {
            return None;
        }
        let out = unsafe { self.data.item_ptr(self.ix).as_mut() };
        self.ix += 1;
        return Some(out);
    }
}

pub struct Iter<'a, T> {
    data: &'a SqVec<T>,
    ix: u32,
}

impl<'a, T> Iterator for Iter<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.ix == self.data.len() as u32 {
            return None;
        }
        let out = unsafe { self.data.item_ptr(self.ix).as_ref() };
        self.ix += 1;
        return Some(out);
    }
}

pub struct TIter<'a, T> {
    seg_ptr: NonNull<T>,
    seg_ix: u32,
    marker: PhantomData<&'a SqVec<T>>,
}

impl<T> Index<u32> for SqVec<T> {
    type Output = T;

    fn index(&self, ix: u32) -> &Self::Output {
        assert!(ix < self.len, "Index out of bounds.");
        unsafe { self.item_ptr(ix).as_ref() }
    }
}

impl<T> IndexMut<u32> for SqVec<T> {
    fn index_mut(&mut self, ix: u32) -> &mut Self::Output {
        assert!(ix < self.len, "Index out of bounds.");
        unsafe { self.item_ptr(ix).as_mut() }
    }
}

impl<T: Debug> Debug for SqVec<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

unsafe impl<T: Send> Send for SqVec<T> {}
unsafe impl<T: Sync> Sync for SqVec<T> {}

impl<T: PartialEq> PartialEq for SqVec<T> {
    fn eq(&self, other: &Self) -> bool {
        if self.len() != other.len() {
            return false;
        }
        self.iter().eq(other.iter())
    }
}

impl<T: PartialEq> PartialEq<Vec<T>> for SqVec<T> {
    fn eq(&self, other: &Vec<T>) -> bool {
        if self.len() != other.len() {
            return false;
        }
        self.iter().eq(other.iter())
    }
}

impl<T: PartialEq> PartialEq<SqVec<T>> for Vec<T> {
    fn eq(&self, other: &SqVec<T>) -> bool {
        other == self
    }
}

#[cfg(test)]
mod tests {
    use super::SqVec;
    use rand::{Rng, random, rng, thread_rng};

    #[test]
    fn single_box_drop() {
        let mut sqvec = SqVec::<Box<u32>>::new();
        sqvec.push(Box::new(32));
        sqvec.pop();
    }

    #[test]
    fn drop_len_lt_cap() {
        let mut sqvec = SqVec::<Box<u32>>::new();
        for i in 0..200 {
            sqvec.push(Box::new(i));
        }
        for _ in 0..170 {
            sqvec.pop();
        }
    }

    #[test]
    fn drop_empty() {
        let _sqvec = SqVec::<Box<u32>>::new();
    }

    #[test]
    /// randomised test of mutating api
    fn fuzz() {
        let mut rng = thread_rng();
        let mut vec = Vec::<Box<u32>>::new();
        let mut sqvec = SqVec::<Box<u32>>::new();
        for _ in 0u32..2048 {
            let choice: u8 = rng.gen_range(0..=4);
            match choice {
                0 => {
                    if !vec.is_empty() {
                        let ix = rng.gen_range(0..vec.len());
                        assert_eq!(vec[ix], sqvec[ix as u32], "Not equal at index: {ix}");
                    }
                }
                1 => {
                    assert_eq!(vec.pop(), sqvec.pop());
                }
                2..=4 => {
                    let val: Box<u32> = Box::new(random());
                    vec.push(val.clone());
                    sqvec.push(val.clone());
                }
                _ => (),
            };
            assert_eq!(vec.len(), sqvec.len() as usize);
        }
    }

    #[test]
    fn lifetime() {
        let mut sqvec = SqVec::<Box<u32>>::new();
        sqvec.push(Box::new(23));
        let a = sqvec.pop().unwrap();
        drop(a);
        sqvec.push(Box::new(324));
    }

    #[test]
    fn a_titer() {
        let mut sqvec = SqVec::<u32>::new();
        for i in 0..2048 {
            sqvec.push(i);
        }
        let mut si = sqvec.iter_t();
        for _ in 0..4 {
            si.next();
            println!("{:?}", &si);
        }
    }

    #[test]
    fn iteration() {
        let mut rng = rng();
        let mut sqvec = SqVec::<u32>::new();
        let mut vec = Vec::<u32>::new();
        for _ in 0u32..2048 {
            let a = rng.random();
            sqvec.push(a);
            vec.push(a);
        }
        let mut sqiter = sqvec.iter();
        let mut sqiter2 = sqvec.iter_t();
        let mut viter = vec.iter();
        for i in 0u32..2048 {
            // assert_eq!(sqiter.next(), viter.next(), "{i}");
            assert_eq!(sqiter.next(), sqiter2.next(), "{i}");
        }
    }

    #[test]
    fn equal() {
        let mut sqvec1: SqVec<u32> = SqVec::new();
        let mut sqvec2: SqVec<u32> = SqVec::new();
        let mut vec: Vec<u32> = Vec::new();

        for i in [0, 5, 1, 56, 2] {
            sqvec1.push(i);
            sqvec2.push(i);
            vec.push(i);
        }

        assert_eq!(sqvec1, sqvec2);
        assert_eq!(sqvec1, vec);
        assert_eq!(vec, sqvec1);

        assert_ne!(sqvec!(3, 5, 2), sqvec!(2, 3));
        assert_ne!(sqvec!(33, 1, 34, 1), vec![23, 3, 5]);
        assert_ne!(sqvec!(2, 3, 6), sqvec!(2, 3, 6, 5));
    }
}
