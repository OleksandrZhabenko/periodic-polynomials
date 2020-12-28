-- |
-- Module      :  Data.Periodic
-- Copyright   :  (c) OleksandrZhabenko 2020
-- License     :  MIT
-- Stability   :  Experimental
-- Maintainer  :  olexandr543@yahoo.com
--
-- A library for working with periodic polynomials (very basic functionality). 
-- Provides also simple tools to make a numerical function a periodic (or somewhat similar) one. 
-- 

module Data.Periodic (
  -- * The simplest finite periodic polynomials
  polyG1
  , trigPolyCos
  , trigPolySin
  , trigPoly
  -- * Periodizer functions
  , periodizer
  , concatPeriodizer
  -- * Periodizer applications
  , polyG2
  , polyG3
) where

import qualified Data.Vector as V

-- | A finite trigonometric polynomial of sines. The 'V.Vector' argument is used to produce its coefficients (weights) by applying to each of the element the function 
-- @f:: a -> a@ given as the first argument.
trigPolySin :: (Floating a) => (a -> a) -> V.Vector a -> a -> a
trigPolySin f = polyG1 f (sin)

-- | A finite trigonometric polynomial of cosines. The 'V.Vector' argument is used to produce its coefficients (weights) by applying to each of the element the function 
-- @f:: a -> a@ given as the first argument.
trigPolyCos :: (Floating a) => (a -> a) -> V.Vector a -> a -> a
trigPolyCos f = polyG1 f (cos)

-- | Sum of the sine and cosine finite trigonometric polynomials. Can represent Fourier series (without the first coefficient), but no numerical high accuracy is guaranteed.
trigPoly :: (Floating a) => (a -> a) -> V.Vector a -> (a -> a) -> V.Vector a -> a -> a
trigPoly f v1 g v2 y 
 | V.null v1 = polyG1 g (sin) v2 y
 | V.null v2 = polyG1 f (cos) v1 y
 | otherwise = (V.sum . V.zipWith (*) (V.map g v2) . V.generate (V.length v2) $ (\i -> sin (fromIntegral (i + 1) * y))) + 
    (V.sum . V.zipWith (*) (V.map f v1) . V.generate (V.length v1) $ (\i -> cos (fromIntegral (i + 1) * y)))

-- | Makes a function @f :: a -> b@ periodic with the period given by the third argument and the starting point given by the second argument. 
periodizer :: (RealFrac a) => (a -> b) -> a -> a -> a -> b
periodizer f x0 period x 
 | period /= 0.0 = let delta = (x - x0) / abs period in f (x0 + (abs period) * (delta - (fromIntegral . truncate $ delta)))
 | otherwise = error "Data.Periodic.periodizer: Not defined for the zero period. "

-- | Modified periodizer that tries to concat the pieces of the function so that it can be (generally speaking) continuous. Needs more mathematical studies. 
concatPeriodizer :: (RealFrac a, Num b) => (a -> b) -> a -> a -> a -> b
concatPeriodizer f x0 period x 
 | period /= 0.0 = 
   let delta = (x - x0) / abs period
       deltaI = truncate delta in f (x0 + (abs period) * (delta - (fromIntegral deltaI))) + fromIntegral deltaI * (f (x0 + abs period) - f x0)
 | otherwise = error "Data.Periodic.concatPeriodizer: Not defined for the zero period. "

-- | The first function @f :: a -> a@ is applied to the vector to produce weighted coefficients for the sum and the second one @g :: a -> a@ is used as a basis 
-- function. For the periodic function g the resulting function is also periodic with the same period. Among possible variants there are finite trigonometric 
-- polynomials. See as examples 'trigPolySin' and 'trigPolyCos' functions.
polyG1 :: (Floating a) => (a -> a) -> (a -> a) -> V.Vector a -> a -> a
polyG1 f g v y 
 | V.null v = 0.0
 | otherwise = V.sum . V.zipWith (*) (V.map f v) . V.generate (V.length v) $ (\i -> g (fromIntegral (i + 1) * y))

-- | Instead of simply use the second function in 'polyG1' it applies to it a 'periodizer' with the given arguments.
polyG2 :: (RealFrac a) => (a -> a) -> (a -> a) -> a -> a -> V.Vector a -> a -> a
polyG2 f g period y0 v y 
 | V.null v = 0.0
 | otherwise = V.sum . V.zipWith (*) (V.map f v) . V.generate (V.length v) $ (\i -> periodizer g y0 period (fromIntegral (i + 1) * y))

-- | Instead of simply use the second function in 'polyG1' it applies to it a 'concatPeriodizer' with the given arguments.
polyG3 :: (RealFrac a) => (a -> a) -> (a -> a) -> a -> a -> V.Vector a -> a -> a
polyG3 f g period y0 v y 
 | V.null v = 0.0
 | otherwise = V.sum . V.zipWith (*) (V.map f v) . V.generate (V.length v) $ (\i -> concatPeriodizer g y0 period (fromIntegral (i + 1) * y))
