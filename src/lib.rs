//! For representing real numbers as fractions using integer types.

use std::fmt::{Debug, Display, Formatter};
use std::ops::{Add, Sub, Mul, Div, Neg};

/// ...
pub trait FractionTerm: PartialOrd + Copy + private::Sealed {
	const MAX_TERM_COUNT: usize;
	fn from_f64(num: f64) -> Self;
	fn into_f64(self) -> f64;
	fn checked_add(self, rhs: Self) -> Option<Self>;
	fn checked_sub(self, rhs: Self) -> Option<Self>;
	fn checked_mul(self, rhs: Self) -> Option<Self>;
	fn checked_div(self, rhs: Self) -> Option<Self>;
	fn gcd(self, divisor: Self) -> Self;
}
mod private {
    pub trait Sealed {}
}

macro_rules! impl_fraction_term {
	($type:ty) => {
		impl private::Sealed for $type {}
		impl FractionTerm for $type {
			const MAX_TERM_COUNT: usize = (1 + ((Self::MAX.ilog2() + 1) / 2)) as usize;
			
			fn from_f64(num: f64) -> Self {
				num as Self
			}
			fn into_f64(self) -> f64 {
				self as f64
			}
			
			fn checked_add(self, rhs: Self) -> Option<Self> {
				self.checked_add(rhs)
			}
			fn checked_sub(self, rhs: Self) -> Option<Self> {
				self.checked_sub(rhs)
			}
			fn checked_mul(self, rhs: Self) -> Option<Self> {
				self.checked_mul(rhs)
			}
			fn checked_div(self, rhs: Self) -> Option<Self> {
				self.checked_div(rhs)
			}
			
			fn gcd(mut self, mut divisor: Self) -> Self {
				while self != 0 {
					let temp = self;
					self = divisor % self;
					divisor = temp;
				}
				divisor
			}
		}
	}
}

impl_fraction_term!(u8);
impl_fraction_term!(u16);
impl_fraction_term!(u32);
impl_fraction_term!(u64);
impl_fraction_term!(u128);

/// ...
#[derive(Debug, Clone)]
pub enum Fraction<T: FractionTerm> {
	Pos(T, T),
	Neg(T, T),
}

impl<T: FractionTerm> Fraction<T> {
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
	
	pub fn signum(&self) -> i8 {
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
	
	pub fn to_f64(&self) -> f64 {
		match self {
			Self::Pos(n, d) =>  n.into_f64() / d.into_f64(),
			Self::Neg(n, d) => -n.into_f64() / d.into_f64(),
		}
	}
}

impl<T: FractionTerm> Add for Fraction<T> {
	type Output = Self;
	fn add(self, rhs: Self) -> Self::Output {
		self.checked_add(rhs).unwrap()
	}
}

impl<T: FractionTerm> Sub for Fraction<T> {
	type Output = Self;
	fn sub(self, rhs: Self) -> Self::Output {
		self.checked_sub(rhs).unwrap()
	}
}

impl<T: FractionTerm> Mul for Fraction<T> {
	type Output = Self;
	fn mul(self, rhs: Self) -> Self::Output {
		self.checked_mul(rhs).unwrap()
	}
}

impl<T: FractionTerm> Div for Fraction<T> {
	type Output = Self;
	fn div(self, rhs: Self) -> Self::Output {
		self.checked_div(rhs).unwrap()
	}
}

impl<T: FractionTerm> Neg for Fraction<T> {
	type Output = Self;
	fn neg(self) -> Self::Output {
		match self {
			Self::Pos(n, d) => Self::Neg(n, d),
			Self::Neg(n, d) => Self::Pos(n, d),
		}
	}
}

impl<T: FractionTerm> From<f32> for Fraction<T> {
	fn from(value: f32) -> Self {
		Fraction::from(value as f64)
	}
}

impl<T: FractionTerm> From<f64> for Fraction<T> {
	fn from(mut num: f64) -> Self {
		/*!
			Returns an approximate fraction for the given number, saturating to the type's bounds.
			
			# Algorithm
			
			```text
			f(x) = 1 / (f(x-1) - r(x-1))
			r(x) = round(f(x))
			f(0) = number
			     = r(0) ± 1/(r(1) ± 1/(r(2) ..))
			
			Ex:
				f(0)=1.55, f(1)=2.222, f(2)=4.5, f(3)=2
				Depending on the integer type, this result may be 31/20, 14/9, 3/2, or 2/1.
				
				        1             1         -9   31
				2 - ————————— = 2 - ————— = 2 + —— = —— = 1.55
				          1             2       20   20
				    2 + —————       2 + —
				            1           9
				        4 + —
				            2
			```
			
			# Examples
			
			```
			use fraction::Fraction;
			assert_eq!("1/3",      format!("{}", Fraction::<u8> ::from(0.333_f64)));
			assert_eq!("333/1000", format!("{}", Fraction::<u16>::from(0.333_f64)));
			assert_eq!("-50/1",    format!("{}", Fraction::<u64>::from(-50_f32)));
			assert_eq!("22/7",     format!("{}", Fraction::<u8> ::from(std::f64::consts::PI)));
			```
		*/
		
		if num.is_nan() {
			return Self::Pos(T::from_f64(0.0), T::from_f64(0.0))
		}
		
		let mut term_list: Vec<Self> = Vec::with_capacity(T::MAX_TERM_COUNT);
		
		let mut limit = T::from_f64(f64::INFINITY).into_f64();
		
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
				if let Some(new_frac) = frac.checked_add(term_list[term_index].clone()) {
					frac = new_frac;
					continue
				}
				continue 'simplify
			}
			return frac
		}
	}
}

impl<T: FractionTerm> From<Fraction<T>> for f64 {
	fn from(value: Fraction<T>) -> Self {
		value.to_f64()
	}
}

impl<T: FractionTerm> From<Fraction<T>> for f32 {
	fn from(value: Fraction<T>) -> Self {
		value.to_f64() as f32
	}
}

macro_rules! impl_from_fraction_for_fraction {
	($F:ty, $T:ty) => {
		impl From<Fraction<$F>> for Fraction<$T>{
			fn from(value: Fraction<$F>) -> Self {
				match value {
					Fraction::<$F>::Pos(n, d) => Self::Pos(n.into(), d.into()),
					Fraction::<$F>::Neg(n, d) => Self::Neg(n.into(), d.into()),
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

impl_from_fraction_for_fraction!(u8,  usize);
impl_from_fraction_for_fraction!(u16, usize);

impl<T: FractionTerm + Display> Display for Fraction<T> {
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
	fn accuracy() {
		fn assert<I: FractionTerm + Display>() {
			for n in 1..=100 {
				for d in 1..=100 {
					let n: f64 = n.into();
					let d: f64 = d.into();
					let f = Fraction::<I>::from(n/d);
					// println!("{n}/{d} -> {f}, {}", n/d);
					let v: f64 = f.into();
					assert!(v > (n-1.0)/d);
					assert!(v < (n+1.0)/d);
				}
			}
			assert_eq!(Fraction::<I>::from(f64::INFINITY).to_f64(), I::from_f64(f64::INFINITY).into_f64());
			assert_eq!(Fraction::<I>::from(f64::NEG_INFINITY).to_f64(), -I::from_f64(f64::INFINITY).into_f64());
			assert!(Fraction::<I>::from(f64::NAN).to_f64().is_nan());
		}
		assert::<u8>();
		assert::<u64>();
		assert::<u128>();
	}
	
	#[test]
	fn arithmetic() {
		println!("{}", Fraction::Neg(5_u8, 0) + Fraction::Pos(6_u8, 0));
		assert_eq!("17/15",  format!("{}", Fraction::Pos(1_u8, 3) + Fraction::Pos(4, 5)));
		assert_eq!("-7/15",  format!("{}", Fraction::Pos(1_u16, 3) - Fraction::Pos(4, 5)));
		assert_eq!("8/6",    format!("{}", Fraction::Pos(4_u128, 6) + Fraction::Pos(2, 3)));
		assert_eq!("112/6",  format!("{}", Fraction::Pos(8_u8, 3) * Fraction::Pos(14, 2)));
		assert_eq!("-10/50", format!("{}", Fraction::Neg(1_u8, 10) / Fraction::Pos(5, 10)));
	}
}