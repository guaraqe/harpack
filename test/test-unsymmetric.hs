import Test.Tasty
import Test.Tasty.HUnit

import qualified Data.Vector.Storable.Mutable as MV
import qualified Data.Vector.Storable as V
import Data.Vector.Storable ((!))
import Numeric.LinearAlgebra.HArpack
import Control.Monad
import Foreign.C.Types
import System.IO.Unsafe
import Data.Complex

import Data.List
import Data.Ord
import Data.Maybe
import Data.Function (on)

main = defaultMain tests

tests :: TestTree
tests = testGroup "Tests" [unitTests]

unitTests = testGroup "Unit tests"
  [ testCase "Diagonal Matrix LM"   $ diagMatLMCheck   @?= True
  , testCase "Diagonal Matrix SM"   $ diagMatSMCheck   @?= True
  , testCase "Non-Symmetric Sanity" $ complexEigsCheck @?= True ]

diagMat :: SparseMatrix
diagMat = SparseMatrix 1000
            [((n-1,n-1),fromIntegral n) | n <- [1 .. 1000]]

complexEigs :: SparseMatrix
complexEigs = SparseMatrix 7
  [((0,0),3)
  ,((0,1),-2)
  ,((1,0),4)
  ,((1,1),-1)
  ,((2,2),2)
  ,((3,3),1)
  ,((4,4),0.5)
  ,((5,5),0.25)
  ,((6,6),0.125)]

diagMatLMCheck :: Bool
diagMatLMCheck =
  let RReal ePairs' = eigs diagMat (Symmetric SLM) 6
      ePairs = decreasingR ePairs'
      evals = map fst ePairs
      trueEvals = [1000,999..]
      diff = map abs $ zipWith (-) evals trueEvals
  in all (< 1e-10) diff

decreasingC = sortBy (compare `on` (\x -> - realPart (abs (fst x))))

decreasingR = sortBy (compare `on` (\x -> - abs (fst x)))

diagMatSMCheck :: Bool
diagMatSMCheck =
  let RComplex ePairs = eigs diagMat (NonSymmetric NSM) 6
      evals = map fst ePairs
      trueEvals = map (:+ 0) [1,2..]
      diff = map (realPart . abs) $ zipWith (-) evals trueEvals
  in all (< 1e-10) diff

complexEigsCheck :: Bool
complexEigsCheck =
  let RComplex ePairs' = eigs complexEigs (NonSymmetric NLM) 5
      ePairs = decreasingC ePairs'
      evals = map fst ePairs
      trueEvals = [1 :+ 2, 1 :+ (-2), 2, 1, 0.5]
      diff = map (realPart . abs) $ zipWith (-) evals trueEvals
  in all (< 1e-10) diff
