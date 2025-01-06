use std::alloc::{self, alloc, dealloc, realloc, Layout};
use std::ops::{Index, IndexMut};
use std::ptr::NonNull;
// use std::ptr::Unique;

#[derive(Debug)]
pub struct SqVec<T> {
    dope: NonNull<NonNull<T>>,
    dope_cap: u32,
    len: u32,
    alloc_seg_count: u32,
    log_seglen: u32,
}

impl<T> SqVec<T> {
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
            log_seglen: 1,
        }
    }

    fn grow_dope(&mut self) {
        let (new_cap, new_layout) = if self.dope_cap == 0 {
            (2, Layout::array::<NonNull<T>>(2).unwrap())
        } else {
            self.log_seglen += 1;
            let new_cap = self.dope_cap + Self::segs_of_len(self.log_seglen);
            let new_layout = Layout::array::<NonNull<T>>(new_cap as usize).unwrap();
            (new_cap, new_layout)
        };

        assert!(
            new_layout.size() <= isize::MAX as usize,
                "Allocation too large"
        );

        let handle = if self.dope_cap == 0 {
            unsafe { alloc(new_layout) }
        } else {
            let old_layout = Layout::array::<NonNull<T>>(self.dope_cap as usize).unwrap();
            let old_handle = self.dope.as_ptr() as *mut u8;
            unsafe { realloc(old_handle, old_layout, new_layout.size()) }
        };

        self.dope = NonNull::new(handle as *mut NonNull<T>)
        .unwrap_or_else(|| alloc::handle_alloc_error(new_layout));
        self.dope_cap = new_cap;
    }

    /// Allocates a new segment of length `1 << self.log_seglen` at index `segnum` in dope.
    fn alloc_seg(&mut self, segnum: u32) {
        let layout = Layout::array::<T>(1 << self.log_seglen).unwrap();
        let handle = unsafe { alloc(layout) };
        let ptr =
        NonNull::new(handle as *mut T).unwrap_or_else(|| alloc::handle_alloc_error(layout));
        unsafe { self.dope.add(segnum as usize).write(ptr) };
        self.alloc_seg_count += 1;
    }

    /// Adds `val` to the end of the SqVec, increasing its length by 1.
    pub fn push(&mut self, val: T) {
        let (segnum, offset) = Self::mapping(self.len);
        if (offset == 0) && (segnum == self.alloc_seg_count) {
            if segnum == self.dope_cap {
                self.grow_dope();
            }
            self.alloc_seg(segnum);
        }
        unsafe { self.item_ptr_inner(segnum, offset).write(val) };
        self.len += 1;
    }

    /// Removes the last value on the SqVec and returns it in an Option.
    /// If the SqVec is empty, Returns None
    pub fn pop(&mut self) -> Option<T> {
        if self.len == 0 {
            return None;
        }
        let (segnum, offset) = Self::mapping(self.len - 1);
        if (segnum == Self::first_seg(self.log_seglen)) && (offset == 0) {
            self.log_seglen = 1.max(self.log_seglen - 1);
        }
        self.len -= 1;
        unsafe {
            let ptr = self.item_ptr_inner(segnum, offset);
            let out = ptr.read();
            ptr.drop_in_place();
            return Some(out);
        }
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.len as usize
    }

    #[inline(always)]
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
    #[inline(always)]
    fn first_seg(log_seglen: u32) -> u32 {
        if log_seglen == 1 {
            0
        } else {
            3 * (1 << log_seglen - 2) - 1
        }
    }

    /// Gives the dope index of the last segment of length (1 << log_seglen)
    #[inline(always)]
    fn last_seg(log_seglen: u32) -> u32 {
        3 * (1 << (log_seglen - 1)) - 2
    }

    #[inline(always)]
    fn segs_of_len(log_seglen: u32) -> u32 {
        if log_seglen == 1 {
            2
        } else {
            3 * (1 << (log_seglen - 2))
        }
    }

    /// Returns the pointer to the item in the SqVec at index `ix` (Does not check if ix is < self.len)
    #[inline(always)]
    unsafe fn item_ptr(&self, ix: u32) -> NonNull<T> {
        let (segnum, offset) = Self::mapping(ix);
        unsafe {self.item_ptr_inner(segnum, offset)}
    }

    #[inline(always)]
    unsafe fn item_ptr_inner(&self, segnum: u32, offset: u32) -> NonNull<T> {
        unsafe {self.dope.add(segnum as usize).read().add(offset as usize)}
    }
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

pub struct Iter<'a, T> {
    ix: u32,
    inner: dyn std::iter::Iterator<Item = &'a T>,
}

impl<'a, T> Iterator for Iter<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        todo!();
    }
}

impl<T> Drop for SqVec<T> {
    fn drop(&mut self) {
        if self.dope_cap == 0 {
            return;
        }
        // drop items in SqVec
        for i in 0..self.len {
            unsafe { self.item_ptr(i).drop_in_place() };
        }

        let dope_slice = unsafe {
            std::slice::from_raw_parts(self.dope.as_ptr(), self.alloc_seg_count as usize)
        };
        let mut ix1 = 0_usize;
        for lsn in 1..=self.log_seglen {
            let ix2 = Self::last_seg(lsn).min(self.alloc_seg_count.saturating_sub(1)) as usize;
            let layout = Layout::array::<T>(1 << lsn).unwrap();
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

#[cfg(test)]
mod tests {
    use super::SqVec;
    use rand::{random, thread_rng, Rng};

    #[test]
    fn fuzz() {
        let mut rng = thread_rng();
        let mut vec = Vec::<u32>::new();
        let mut sqvec = SqVec::<u32>::new();
        for _ in 0..(1 << 16) {
            let choice: u8 = rng.gen_range(0..=3);
            match choice {
                0 => {
                    if !vec.is_empty() {
                        let ix = rng.gen_range(0..vec.len());
                        assert_eq!(vec[ix], sqvec[ix as u32]);
                    }
                }
                1 => {
                    assert_eq!(vec.pop(), sqvec.pop());
                }
                2..=3 => {
                    let val = random();
                    vec.push(val);
                    sqvec.push(val);
                }
                _ => (),
            };
            assert_eq!(vec.len(), sqvec.len());
        }
    }
}
