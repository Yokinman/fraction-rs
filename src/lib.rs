//! Representing real numbers as fractions using integer types.

use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::ops::{Add, Sub, Mul, Div, Neg};

/// ...
pub trait FracTerm: PartialOrd + Copy + private::Sealed {
	const MAX_TERM_COUNT: usize;
	
	fn from_f64(num: f64) -> Self;
	fn to_f64(self) -> f64;
	
	fn checked_add(self, rhs: Self) -> Option<Self>;
	fn checked_sub(self, rhs: Self) -> Option<Self>;
	fn checked_mul(self, rhs: Self) -> Option<Self>;
	fn checked_div(self, rhs: Self) -> Option<Self>;
	fn checked_rem(self, rhs: Self) -> Option<Self>;
	
	fn gcd(mut self, mut other: Self) -> Self {
		while let Some(rem) = self.checked_rem(other) {
			self = other;
			other = rem;
		}
		self
	}
}
mod private {
    pub trait Sealed {}
}

macro_rules! impl_fraction_term {
	($type:ty) => {
		impl private::Sealed for $type {}
		impl FracTerm for $type {
			const MAX_TERM_COUNT: usize = (1 + ((Self::MAX.ilog2() + 1) / 2)) as usize;
			
			fn from_f64(num: f64) -> Self { num as Self }
			fn to_f64(self) -> f64 { self as f64 }
			
			fn checked_add(self, rhs: Self) -> Option<Self> { self.checked_add(rhs) }
			fn checked_sub(self, rhs: Self) -> Option<Self> { self.checked_sub(rhs) }
			fn checked_mul(self, rhs: Self) -> Option<Self> { self.checked_mul(rhs) }
			fn checked_div(self, rhs: Self) -> Option<Self> { self.checked_div(rhs) }
			fn checked_rem(self, rhs: Self) -> Option<Self> { self.checked_rem(rhs) }
		}
	}
}

impl_fraction_term!(u8);
impl_fraction_term!(u16);
impl_fraction_term!(u32);
impl_fraction_term!(u64);
impl_fraction_term!(u128);

/// ...
#[derive(Copy, Clone, Debug)]
pub enum Frac<T: FracTerm> {
	Pos(T, T),
	Neg(T, T),
}

impl<T: FracTerm> Frac<T> {
	pub fn numer(&self) -> T {
		match *self {
			Self::Pos(n, _) => n,
			Self::Neg(n, _) => n,
		}
	}
	
	pub fn denom(&self) -> T {
		match *self {
			Self::Pos(_, d) => d,
			Self::Neg(_, d) => d,
		}
	}
	
	pub fn signum(&self) -> i32 {
		match self {
			Self::Pos(..) =>  1,
			Self::Neg(..) => -1,
		}
	}
	
	pub fn recip(mut self) -> Self {
		match &mut self {
			Self::Pos(n, d) => std::mem::swap(n, d),
			Self::Neg(n, d) => std::mem::swap(n, d),
		}
		self
	}
	
	pub fn checked_add(self, frac: Self) -> Option<Self> {
		let d_gcd = self.denom().gcd(frac.denom());
		let rhs_m = self.denom().checked_div(d_gcd).unwrap_or(T::from_f64(1.0));
		let rhs_n = frac.numer().checked_mul(rhs_m)?;
		let lhs_m = frac.denom().checked_div(d_gcd).unwrap_or(T::from_f64(1.0));
		let lhs_n = self.numer().checked_mul(lhs_m)?;
		let new_d = self.denom().checked_mul(lhs_m)?;
		Some(
			match (self, frac) {
				(Self::Pos(..), Self::Pos(..)) => Self::Pos(lhs_n.checked_add(rhs_n)?, new_d),
				(Self::Neg(..), Self::Neg(..)) => Self::Neg(lhs_n.checked_add(rhs_n)?, new_d),
				(Self::Pos(..), Self::Neg(..)) => {
					if lhs_n >= rhs_n {
						Self::Pos(lhs_n.checked_sub(rhs_n)?, new_d)
					} else {
						Self::Neg(rhs_n.checked_sub(lhs_n)?, new_d)
					}
				},
				(Self::Neg(..), Self::Pos(..)) => {
					if lhs_n >= rhs_n {
						Self::Neg(lhs_n.checked_sub(rhs_n)?, new_d)
					} else {
						Self::Pos(rhs_n.checked_sub(lhs_n)?, new_d)
					}
				},
			}
		)
	}
	
	pub fn checked_sub(self, frac: Self) -> Option<Self> {
		self.checked_add(-frac)
	}
	
	pub fn checked_mul(self, frac: Self) -> Option<Self> {
		let n = self.numer().checked_mul(frac.numer())?;
		let d = self.denom().checked_mul(frac.denom())?;
		if self.signum() == frac.signum() {
			Some(Self::Pos(n, d))
		} else {
			Some(Self::Neg(n, d))
		}
	}
	
	pub fn checked_div(self, frac: Self) -> Option<Self> {
		self.checked_mul(frac.recip())
	}
}

impl<T: FracTerm> Default for Frac<T> {
	fn default() -> Self {
		Self::Pos(T::from_f64(1.0), T::from_f64(1.0))
	}
}

impl<T: FracTerm> PartialEq<Self> for Frac<T> {
	fn eq(&self, other: &Self) -> bool {
		self.partial_cmp(other) == Some(Ordering::Equal)
	}
}

impl<T: FracTerm + PartialOrd> PartialOrd<Self> for Frac<T> {
	//! - 0/0 != 0/0
	//! - A/N <> B/N === A <> B
	//! - A/B <> C/D === AD <> CB
	
	fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
		let (mut s_n, mut s_d) = self.into();
		let (mut o_n, mut o_d) = other.into();
		
		 // Zero Comparison (±0/0 != ±0/0, ±0/A == ±0/B)
		let zero = T::from_f64(0.0);
		if s_n == zero && o_n == zero {
			return if s_n == s_d || o_n == o_d {
				None
			} else {
				Some(Ordering::Equal)
			}
		}
		
		 // Sign Comparison (-A/B <= +C/D):
		let sign_order = self.signum().cmp(&other.signum());
		if sign_order != Ordering::Equal {
			return Some(sign_order)
		}
		
		 // Negative Swap (-A/B <> -C/D === A/B <> C/D):
		if let Self::Neg(..) = self {
			std::mem::swap(&mut s_n, &mut o_n);
			std::mem::swap(&mut s_d, &mut o_d);
		}
		
		 // Numerator Comparison (A/N <> B/N === A <> B):
		if s_d == o_d {
			return s_n.partial_cmp(&o_n)
		}
		
		 // Continued Comparison:
		loop {
			// 22/7 <> 18/5 => 3[1/7] <> 3[3/5] (3=3, try remainder reciprocals)
			//  5/3 <>  7/1 => 1[1/3] <> 7[0/1] (1<7 === 22/7 < 18/5)
			let s_r = s_n.checked_div(s_d);
			let o_r = o_n.checked_div(o_d);
			if let (Some(s_r), Some(o_r)) = (s_r, o_r) {
				if s_r != o_r {
					break s_r.partial_cmp(&o_r)
				}
			} else {
				break o_r.partial_cmp(&s_r)
			}
			s_n = std::mem::replace(&mut o_d, s_n.checked_rem(s_d)?);
			o_n = std::mem::replace(&mut s_d, o_n.checked_rem(s_n)?);
		}
	}
}

impl<T: FracTerm> Add for Frac<T> {
	type Output = Self;
	fn add(self, rhs: Self) -> Self::Output {
		self.checked_add(rhs)
			.expect("overflow when adding fractions")
	}
}

impl<T: FracTerm> Sub for Frac<T> {
	type Output = Self;
	fn sub(self, rhs: Self) -> Self::Output {
		self.checked_sub(rhs)
			.expect("overflow when subtracting fractions")
	}
}

impl<T: FracTerm> Mul for Frac<T> {
	type Output = Self;
	fn mul(self, rhs: Self) -> Self::Output {
		self.checked_mul(rhs)
			.expect("overflow when multiplying fractions")
	}
}

impl<T: FracTerm> Div for Frac<T> {
	type Output = Self;
	fn div(self, rhs: Self) -> Self::Output {
		self.checked_div(rhs)
			.expect("overflow when dividing fractions")
	}
}

impl<T: FracTerm> Neg for Frac<T> {
	type Output = Self;
	fn neg(self) -> Self::Output {
		match self {
			Self::Pos(n, d) => Self::Neg(n, d),
			Self::Neg(n, d) => Self::Pos(n, d),
		}
	}
}

impl<T: FracTerm> Mul<i32> for Frac<T> {
	type Output = Self;
	fn mul(mut self, rhs: i32) -> Self::Output {
		let (n, _) = <(&mut T, _)>::from(&mut self);
		
		*n = n.checked_mul(T::from_f64(rhs.abs() as f64))
			.expect("overflow when multiplying fraction by scalar");
		
		if rhs < 0 {
			-self
		} else {
			self
		}
	}
} 

impl<T: FracTerm> Mul<Frac<T>> for i32 {
	type Output = Frac<T>;
	fn mul(self, rhs: Frac<T>) -> Self::Output {
		rhs * self
	}
}

impl<T: FracTerm> From<f32> for Frac<T> {
	//! Generates the closest fraction within the integer type's bounds.
	
	fn from(value: f32) -> Self {
		Frac::from(value as f64)
	}
}

impl<T: FracTerm> From<f64> for Frac<T> {
	//! Generates the closest fraction within the integer type's bounds.
	//! 
	//! # Examples
	//! 
	//! ```
	//! use fraction::Frac;
	//! assert_eq!("1/3",      format!("{}", Frac::<u8> ::from(0.333_f64)));
	//! assert_eq!("333/1000", format!("{}", Frac::<u16>::from(0.333_f64)));
	//! assert_eq!("-50/1",    format!("{}", Frac::<u64>::from(-50_f32)));
	//! assert_eq!("22/7",     format!("{}", Frac::<u8> ::from(std::f64::consts::PI)));
	//! ```
	
	fn from(mut num: f64) -> Self {
		// # Algorithm
		// 
		// f(x) = 1 / (f(x-1) - r(x-1))
		// r(x) = round(f(x))
		// f(0) = number
		//      = r(0) ± 1/(r(1) ± 1/(r(2) ..))
		// 
		// Ex:
		// 
		// f(0) = 4.0155, f(1) = 66.66.., f(2) = -3
		// In this case, the result is either 803/200 <u16> or 4/1 <u8>.
		// 
		//        1          3    803
		// 4 + —————— = 4 + ——— = ——— = 4.015
		//          1       200   200
		//     67 - —
		//          3
		
		if num.is_nan() {
			return Self::Pos(T::from_f64(0.0), T::from_f64(0.0))
		}
		
		let mut term_list: Vec<Self> = Vec::with_capacity(T::MAX_TERM_COUNT);
		
		let mut limit = T::from_f64(f64::INFINITY).to_f64();
		
		let one = T::from_f64(1.0);
		
		'generate: for _ in 0..T::MAX_TERM_COUNT {
			let int = num.round();
			if int.is_sign_positive() {
				term_list.push(Self::Pos(T::from_f64(int), one));
				limit /= int;
			} else {
				term_list.push(Self::Neg(T::from_f64(-int), one));
				limit /= -int;
			}
			if num == int || limit <= 1.0 {
				break 'generate
			}
			num = (num - int).recip();
		}
		
		'simplify: loop {
			let mut frac = term_list.pop().unwrap();
			let mut term_index = term_list.len();
			while term_index != 0 {
				term_index -= 1;
				frac = frac.recip();
				if let Some(new_frac) = frac.checked_add(term_list[term_index]) {
					frac = new_frac;
					continue
				}
				continue 'simplify
			}
			return frac
		}
	}
}

impl<T: FracTerm> From<&Frac<T>> for f64 {
	fn from(value: &Frac<T>) -> Self {
		match *value {
			Frac::Pos(n, d) =>  n.to_f64() / d.to_f64(),
			Frac::Neg(n, d) => -n.to_f64() / d.to_f64(),
		}
	}
}

impl<T: FracTerm> From<Frac<T>> for f32 {
	fn from(value: Frac<T>) -> Self {
		<f64>::from(&value) as f32
	}
}

impl<T: FracTerm> From<Frac<T>> for f64 {
	fn from(value: Frac<T>) -> Self {
		<f64>::from(&value)
	}
}

impl<'a, T: FracTerm> From<&'a mut Frac<T>> for (&'a mut T, &'a mut T) {
	fn from(value: &'a mut Frac<T>) -> Self {
		match *value {
			Frac::Pos(ref mut n, ref mut d) => (n, d),
			Frac::Neg(ref mut n, ref mut d) => (n, d),
		}
	}
}

impl<T: FracTerm> From<Frac<T>> for (T, T) {
	fn from(value: Frac<T>) -> Self {
		match value {
			Frac::Pos(n, d) => (n, d),
			Frac::Neg(n, d) => (n, d),
		}
	}
}

impl<T: FracTerm> From<&Frac<T>> for (T, T) {
	fn from(value: &Frac<T>) -> Self {
		match *value {
			Frac::Pos(n, d) => (n, d),
			Frac::Neg(n, d) => (n, d),
		}
	}
}

impl<T: FracTerm> From<&mut Frac<T>> for (T, T) {
	fn from(value: &mut Frac<T>) -> Self {
		match *value {
			Frac::Pos(n, d) => (n, d),
			Frac::Neg(n, d) => (n, d),
		}
	}
}

macro_rules! impl_from_fraction_for_fraction {
	($F:ty, $T:ty) => {
		impl From<Frac<$F>> for Frac<$T>{
			fn from(value: Frac<$F>) -> Self {
				match value {
					Frac::<$F>::Pos(n, d) => Self::Pos(n.into(), d.into()),
					Frac::<$F>::Neg(n, d) => Self::Neg(n.into(), d.into()),
				}
			}
		}
	}
}

impl_from_fraction_for_fraction!(u8,  u16);

impl_from_fraction_for_fraction!(u8,  u32);
impl_from_fraction_for_fraction!(u16, u32);

impl_from_fraction_for_fraction!(u8,  u64);
impl_from_fraction_for_fraction!(u16, u64);
impl_from_fraction_for_fraction!(u32, u64);

impl_from_fraction_for_fraction!(u8,  u128);
impl_from_fraction_for_fraction!(u16, u128);
impl_from_fraction_for_fraction!(u32, u128);
impl_from_fraction_for_fraction!(u64, u128);

impl<T: FracTerm + Display> Display for Frac<T> {
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		match self {
			Self::Pos(n, d) => write!(f,  "{}/{}", n, d),
			Self::Neg(n, d) => write!(f, "-{}/{}", n, d),
		}
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	
	#[test]
	fn conversion() {
		fn assert<I: FracTerm + Display>() {
			for n in 1..=100 {
				for d in 1..=100 {
					let n: f64 = n.into();
					let d: f64 = d.into();
					let f = Frac::<I>::from(n/d);
					// println!("{n}/{d} -> {f}, {}", n/d);
					let v: f64 = f.into();
					assert!(v > (n-1.0)/d);
					assert!(v < (n+1.0)/d);
				}
			}
			assert_eq!(<f64>::from(Frac::<I>::from(f64::INFINITY)), I::from_f64(f64::INFINITY).to_f64());
			assert_eq!(<f64>::from(Frac::<I>::from(f64::NEG_INFINITY)), -I::from_f64(f64::INFINITY).to_f64());
			assert!(<f64>::from(Frac::<I>::from(f64::NAN)).is_nan());
		}
		assert::<u8>();
		assert::<u64>();
		assert::<u128>();
	}
	
	#[test]
	fn arithmetic() {
		assert_eq!("17/15",  format!("{}", Frac::Pos(1_u8, 3) + Frac::Pos(4, 5)));
		assert_eq!("-7/15",  format!("{}", Frac::Pos(1_u16, 3) - Frac::Pos(4, 5)));
		assert_eq!("8/6",    format!("{}", Frac::Pos(4_u128, 6) + Frac::Pos(2, 3)));
		assert_eq!("112/6",  format!("{}", Frac::Pos(8_u8, 3) * Frac::Pos(14, 2)));
		assert_eq!("-10/50", format!("{}", Frac::Neg(1_u8, 10) / Frac::Pos(5, 10)));
		assert_eq!("2/0",    format!("{}", Frac::Neg(4_u8, 0) + Frac::Pos(6_u8, 0)));
		assert_eq!("-34/15", format!("{}", -2 * Frac::Pos(17_u8, 15)));
		assert_eq!("-0/15",  format!("{}", Frac::Neg(17_u8, 15) * 0));
	}
	
	#[test]
	fn comparison() {
		assert_eq!(Frac::Pos(5_u8, 5), Frac::Pos(5_u8, 5));
		assert_eq!(Frac::Neg(1_u8, 10), Frac::Neg(2_u8, 20));
		assert!(Frac::Pos(1_u8, 10) > Frac::Neg(2_u8, 20));
		assert!(Frac::Neg(1_u8, 10) < Frac::Pos(2_u8, 20));
		assert!(Frac::Neg(1_u8, 10) > Frac::Neg(3_u8, 20));
		assert!(Frac::Neg(1_u8, 10) < Frac::Neg(1_u8, 20));
		assert_eq!(Frac::Pos(u128::MAX - 1, u128::MAX), Frac::Pos(u128::MAX - 1, u128::MAX));
		assert!(Frac::Pos(u128::MAX - 1, u128::MAX) < Frac::Pos(u128::MAX - 1, u128::MAX - 1));
		assert_eq!(Frac::Pos(1_u8, 0), Frac::Pos(1_u8, 0));
		assert!(Frac::Pos(2_u8, 0) > Frac::Pos(1_u8, 0));
		assert!(Frac::Pos(1_u8, 0) > Frac::Pos(1_u8, 1));
		assert_eq!(Frac::Pos(0_u8, 1), Frac::Pos(0_u8, 2));
		assert_eq!(Frac::Pos(0_u8, 1), Frac::Neg(0_u8, 1));
		assert_ne!(Frac::Pos(0_u8, 0), Frac::Pos(0_u8, 0));
		assert_ne!(Frac::Pos(0_u8, 0), Frac::Neg(0_u8, 0));
	}
}