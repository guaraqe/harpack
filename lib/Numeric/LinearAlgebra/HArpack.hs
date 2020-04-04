{-# LANGUAGE ForeignFunctionInterface #-}

module Numeric.LinearAlgebra.HArpack
       ( eigs
       , eigsN
       , eigsS
       , linearOperator
       , SparseMatrix (..)
       , Problem (..)
       , SWhich (..)
       , NWhich (..)
       , Result (..)
       , Dimension
       , HowMany
       , Tolerance
       , Iterations ) where

import           Control.Monad
import           Control.Monad.Loops
import           Data.Complex
import           Data.Vector.Storable         (Vector)
import qualified Data.Vector.Storable         as V
import           Data.Vector.Storable.Mutable (IOVector)
import qualified Data.Vector.Storable.Mutable as MV
import           Foreign
import           Foreign.C.String
import           Foreign.C.Types
import           Foreign.ForeignPtr.Unsafe    as FPUnsafe
import           System.IO.Unsafe

-- | Functions of type 'LinearOperator' represent a matrix vector product.
-- The values in the second vector should be overwritten by operator applied to
-- the first.

data SparseMatrix =
  SparseMatrix { dim ::Int
               , indexes :: [((Int,Int),Double)] }
  deriving Show

data LinearOperator =
  LinearOperator
    { applyOperator :: IOVector CDouble -> IOVector CDouble -> IO () }

linearOperator :: [((Int,Int),Double)] -> LinearOperator
linearOperator l = LinearOperator $ \x y ->
  do MV.set y 0
     forM_ l $ \((i,j),v) -> do
       v' <- (* realToFrac v) <$> MV.read x j
       MV.modify y (+ v') i

--------------------------------------------------------------------------------
-- | Which eigenvalues to calculate:
-- - The first character refers to 'S'ymmetric or 'N'on-Symmetric;
-- - The second refers to 'L'argest or 'S'mallest;
-- - The third refers to 'M'agnitude, 'A'lgebraic, 'R'eal, 'I'mmaginary;
-- - The option 'BE' means 'BE'tween, for central eigenvalues.

data SWhich = SLM | SSM | SLA | SSA | SBE

instance Show SWhich where
  show SLM = "LM"
  show SSM = "SM"
  show SLA = "LA"
  show SSA = "SA"
  show SBE = "BE"

----------

data NWhich = NLM | NSM | NLR | NSR | NLI | NSI

instance Show NWhich where
  show NLM = "LM"
  show NSM = "SM"
  show NLR = "LR"
  show NSR = "SR"
  show NLI = "LI"
  show NSI = "SI"

----------

data Problem = Symmetric SWhich | NonSymmetric NWhich
     deriving (Show)

data Eigenpairs = Error String
                | EReal [(Double,Vector Double)]
                | EComplex [(Complex Double, Vector (Complex Double))]

--------------------------------------------------------------------------------

-- | Encapsulates the arrays and data needed to invoke _naupd routines in ARPACK

data ArpackSetup =
  ArpackSetup { aIDO    :: Ptr CInt
              , aBMAT   :: Ptr CChar
              , aN      :: Ptr CInt
              , aWHICH  :: Ptr CChar
              , aNEV    :: Ptr CInt
              , aTOL    :: Ptr CDouble
              , aRESID  :: Ptr CDouble
              , aNCV    :: Ptr CInt
              , aV      :: Ptr CDouble
              , aLDV    :: Ptr CInt
              , aIPARAM :: Ptr CInt
              , aIPNTR  :: Ptr CInt
              , aWORKD  :: Ptr CDouble
              , aWORKL  :: Ptr CDouble
              , aLWORKL :: Ptr CInt
              , aINFO   :: Ptr CInt }
  deriving (Show)

type Dimension = Int
type HowMany = Int
type Tolerance = Double
type Iterations = Int

-- | The 'setupArpack' function allocates the arrays and data
-- needed by ARPACK to solve an eigenvalue problem

setupArpack :: Problem
            -> Dimension
            -> HowMany
            -> Tolerance
            -> Iterations
            -> IO ArpackSetup

setupArpack (NonSymmetric which') n nev tol maxIter =
  do let ncv = min n (2 * nev + 1)
         lWorkl = 3 * ncv * ncv + 6 * ncv
         ldv = n

     idoPtr <- ptrInt 0
     nPtr   <- ptrInt n
     nevPtr <- ptrInt nev
     ncvPtr <- ptrInt ncv
     ldvPtr <- ptrInt ldv

     bmatPtr  <- ptrCharArray "I"
     whichPtr <- ptrCharArray (show which')

     tolPtr   <- ptrDouble tol

     residPtr <- mallocDoubleArray n
     vPtr     <- mallocDoubleArray (n * ncv)
     workdPtr <- mallocDoubleArray (3 * n)
     worklPtr <- mallocDoubleArray lWorkl

     iParamPtr <- mallocIntArray 11

     pokeElemOff iParamPtr 0 1
     pokeElemOff iParamPtr 2 (fromIntegral maxIter)
     pokeElemOff iParamPtr 3 1
     pokeElemOff iParamPtr 6 1

     iPntrPtr  <- mallocIntArray 14

     lworklPtr <- ptrInt lWorkl
     infoPtr   <- ptrInt 0

     setup <- return $
              ArpackSetup idoPtr bmatPtr nPtr whichPtr nevPtr tolPtr
                          residPtr ncvPtr vPtr ldvPtr iParamPtr iPntrPtr
                          workdPtr worklPtr lworklPtr infoPtr

     dnaupd setup
     return setup

setupArpack (Symmetric which) n nev tol maxIter =
  do let ncv = min n (2 * nev + 1)
         lWorkl = ncv * ncv + 8 * ncv
         ldv = n

     idoPtr <- ptrInt 0
     nPtr   <- ptrInt n
     nevPtr <- ptrInt nev
     ncvPtr <- ptrInt ncv
     ldvPtr <- ptrInt ldv

     bmatPtr  <- ptrCharArray "I"
     whichPtr <- ptrCharArray (show which)

     tolPtr   <- ptrDouble tol

     residPtr <- mallocDoubleArray n
     vPtr     <- mallocDoubleArray (n * ncv)
     workdPtr <- mallocDoubleArray (3 * n)
     worklPtr <- mallocDoubleArray lWorkl

     iParamPtr <- mallocIntArray 11

     pokeElemOff iParamPtr 0 1
     pokeElemOff iParamPtr 2 (fromIntegral maxIter)
     pokeElemOff iParamPtr 3 1
     pokeElemOff iParamPtr 6 1

     iPntrPtr  <- mallocIntArray 11

     lworklPtr <- ptrInt lWorkl
     infoPtr   <- ptrInt 0

     setup <- return $
              ArpackSetup idoPtr bmatPtr nPtr whichPtr nevPtr tolPtr
                          residPtr ncvPtr vPtr ldvPtr iParamPtr iPntrPtr
                          workdPtr worklPtr lworklPtr infoPtr

     dsaupd setup

     return setup

-- |Invoke ARPACK's iteration through the reverse communication interface.
iterateArpack :: Problem -> LinearOperator -> ArpackSetup -> IO ()
iterateArpack (NonSymmetric _) f setup = do
  workdXElem <- peekElemOff (aIPNTR setup) 0
  workdYElem <- peekElemOff (aIPNTR setup) 1

  n <- peek (aN setup)

  xPtr <- newForeignPtr_ $ advancePtr (aWORKD setup) (fromIntegral workdXElem - 1)
  yPtr <- newForeignPtr_ $ advancePtr (aWORKD setup) (fromIntegral workdYElem - 1)

  let x = MV.unsafeFromForeignPtr xPtr 0 (fromIntegral n)
      y = MV.unsafeFromForeignPtr yPtr 0 (fromIntegral n)

  applyOperator f x y

  dnaupd setup

iterateArpack (Symmetric _) f setup = do
  workdXElem <- peekElemOff (aIPNTR setup) 0
  workdYElem <- peekElemOff (aIPNTR setup) 1

  n <- peek (aN setup)

  xPtr <- newForeignPtr_ $ advancePtr (aWORKD setup) (fromIntegral workdXElem - 1)
  yPtr <- newForeignPtr_ $ advancePtr (aWORKD setup) (fromIntegral workdYElem - 1)

  let x = MV.unsafeFromForeignPtr xPtr 0 (fromIntegral n)
      y = MV.unsafeFromForeignPtr yPtr 0 (fromIntegral n)

  applyOperator f x y

  dsaupd setup


arpack :: LinearOperator -> Problem -> Dimension -> HowMany -> Tolerance ->
          Iterations -> IO ArpackSetup
arpack op problem n nev tol mxItr = do
  ar <- setupArpack problem n nev tol mxItr
  iterateArpack problem op ar
  whileM_ ((== 1) <$> peek (aIDO ar)) (iterateArpack problem op ar)
  return ar

--------------------------------------------------------------------------------

data D = DReal (IOVector CDouble)
       | DComplex (IOVector CDouble) (IOVector CDouble)

data Sigma = SReal (Ptr CDouble)
           | SComplex (Ptr CDouble) (Ptr CDouble)

data ArpackResults =
  ArpackResults { eRVEC    :: Ptr CInt
                , eHOWMNY  :: Ptr CChar
                , eSELECT  :: Ptr CInt
                , eD       :: D
                , eZ       :: IOVector CDouble
                , eLDZ     :: Ptr CInt
                , eSIGMA   :: Sigma
                , eWORKEV  :: Maybe (Ptr CDouble)
                , eSetup   :: ArpackSetup }

errCheck :: ArpackResults -> IO Bool
errCheck ar = do
  info <- (peek . aINFO . eSetup) ar
  print info
  return (info < 0)

parseArpackOutput :: Problem -> ArpackSetup -> IO ArpackResults
parseArpackOutput (NonSymmetric _) setup = do
  ncv <- peek (aNCV setup)
  nev <- peek (aNEV setup)
  n <- peek (aN setup)

  rvecPtr <- ptrInt 1
  howmnyPtr <- ptrCharArray "A"
  drVec <- newVectorDouble (nev + 1)
  diVec <- newVectorDouble (nev + 1)
  selectPtr <- mallocIntArray (fromIntegral ncv)
  zVec <- newVectorDouble (n*(nev+1))
  ldzPtr <- ptrInt (fromIntegral n)
  sigmarPtr <- mallocDouble
  sigmaiPtr <- mallocDouble
  workevPtr <- mallocDoubleArray (fromIntegral (3*ncv))

  let (drfPtr,_,_) = MV.unsafeToForeignPtr drVec
      drPtr = FPUnsafe.unsafeForeignPtrToPtr drfPtr
      (difPtr,_,_) = MV.unsafeToForeignPtr diVec
      diPtr = FPUnsafe.unsafeForeignPtrToPtr difPtr
      (zfPtr,_,_) = MV.unsafeToForeignPtr zVec
      zPtr = FPUnsafe.unsafeForeignPtrToPtr zfPtr

  c_dneupd rvecPtr howmnyPtr selectPtr drPtr diPtr zPtr ldzPtr sigmarPtr
    sigmaiPtr workevPtr (aBMAT setup) (aN setup) (aWHICH setup) (aNEV setup)
    (aTOL setup) (aRESID setup) (aNCV setup) (aV setup) (aLDV setup)
    (aIPARAM setup) (aIPNTR setup) (aWORKD setup) (aWORKL setup) (aLWORKL setup)
    (aINFO setup)

  return $ ArpackResults rvecPtr howmnyPtr selectPtr (DComplex drVec diVec) zVec
               ldzPtr (SComplex sigmarPtr sigmaiPtr) (Just workevPtr) setup

parseArpackOutput (Symmetric _) setup = do
  ncv <- peek (aNCV setup)
  nev <- peek (aNEV setup)
  n <- peek (aN setup)

  rvecPtr <- ptrInt 1
  howmnyPtr <- ptrCharArray "A"
  dVec <- newVectorDouble nev
  selectPtr <- mallocIntArray (fromIntegral ncv)
  zVec <- newVectorDouble (n*(nev+1))
  ldzPtr <- ptrInt (fromIntegral n)
  sigmaPtr <- mallocDouble

  let (drfPtr,_,_) = MV.unsafeToForeignPtr dVec
      dPtr = FPUnsafe.unsafeForeignPtrToPtr drfPtr
      (zfPtr,_,_) = MV.unsafeToForeignPtr zVec
      zPtr = FPUnsafe.unsafeForeignPtrToPtr zfPtr

  c_dseupd rvecPtr howmnyPtr selectPtr dPtr zPtr ldzPtr
    sigmaPtr (aBMAT setup) (aN setup) (aWHICH setup) (aNEV setup)
    (aTOL setup) (aRESID setup) (aNCV setup) (aV setup) (aLDV setup)
    (aIPARAM setup) (aIPNTR setup) (aWORKD setup) (aWORKL setup) (aLWORKL setup)
    (aINFO setup)

  return $ ArpackResults rvecPtr howmnyPtr selectPtr (DReal dVec) zVec
               ldzPtr (SReal sigmaPtr) Nothing setup

getComplexEigenvalue :: ArpackResults -> Int -> IO (Complex Double)
getComplexEigenvalue results index = do
  let DComplex reEigs imEigs = eD results
  re <- MV.read reEigs index
  im <- MV.read imEigs index
  return (realToFrac re :+ realToFrac im)


getRealEigenvector :: ArpackResults -> Int -> IO (Vector (Complex Double))
getRealEigenvector results index = do
  n <- (peek . aN . eSetup) results
  ldz <- peek $ eLDZ results
  let z = eZ results
  evec <- V.freeze $ MV.slice (index * fromIntegral ldz) (fromIntegral n) z
  return $ V.map ((:+0) . realToFrac) evec

getComplexEigenvector
  :: ArpackResults -> Int -> Bool -> IO (Vector (Complex Double))
getComplexEigenvector results index conjugate = do
  n <- (peek . aN . eSetup) results
  ldz <- peek $ eLDZ results
  let z = eZ results
  reEvec <- V.freeze $ MV.slice (index * fromIntegral ldz) (fromIntegral n) z
  imEvec <- V.freeze $ MV.slice ((index+1) * fromIntegral ldz) (fromIntegral n) z
  if conjugate
     then return $
          V.zipWith (:+) (V.map realToFrac reEvec)
                         (V.map (negate . realToFrac) imEvec)
     else return $
          V.zipWith (:+) (V.map realToFrac reEvec)
                         (V.map realToFrac imEvec)

getEigenpair
  :: ArpackResults -> Int -> IO (Complex Double, V.Vector (Complex Double))
getEigenpair results index = do
  eig <- getComplexEigenvalue results index
  let evec = case compare (imagPart eig) 0 of
        EQ -> getRealEigenvector results index
        LT -> getComplexEigenvector results (index - 1) True
        GT -> getComplexEigenvector results index False
  evec' <- evec
  return (eig,evec')

-- | The function 'eigs' computes a few eigenvalues of a real valued
-- unsymmetric matrix by invoking the ARPACK library. First element is true if
-- ARPACK converged. If converged the second element is a list of eigenpairs.

eigsN :: LinearOperator
      -> NWhich
      -> Dimension
      -> HowMany
      -> Tolerance
      -> Iterations
      -> Maybe [(Complex Double, Vector (Complex Double))]
eigsN f which' n' nev' tol' iters = unsafePerformIO $ do
  setup <- arpack f (NonSymmetric which') n' nev' tol' iters
  results <- parseArpackOutput (NonSymmetric which') setup
  err <- errCheck results
  if err
    then return Nothing
    else do
      pairs <- mapM (getEigenpair results) [0..(nev'-1)]
      return (Just pairs)

data Result = RError
            | RReal [(Double,Vector Double)]
            | RComplex [(Complex Double, Vector (Complex Double))]
              deriving Show

eigs :: SparseMatrix
     -> Problem
     -> HowMany
     -> Result
eigs m (NonSymmetric w) h =
  case eigsN (linearOperator (indexes m)) w (dim m) h (1e-10) 6000 of
  Nothing -> RError
  Just x -> RComplex x
eigs m (Symmetric w) h =
  case eigsS (linearOperator (indexes m)) w (dim m) h (1e-10) 6000 of
  Nothing -> RError
  Just x -> RReal x

--------------------------------------------------------------------------------

getSymmetricEigenvalue :: ArpackResults -> Int -> IO Double
getSymmetricEigenvalue results index = do
  let DReal eigV = eD results
  eig <- MV.read eigV index
  return (realToFrac eig)

getSymmetricEigenvector :: ArpackResults -> Int -> IO (Vector Double)
getSymmetricEigenvector results index = do
  n <- fromIntegral <$> (peek . aN . eSetup) results
  ldz <- peek $ eLDZ results
  let z = eZ results
  evec <- V.freeze $ MV.slice (index * fromIntegral ldz) n z
  return $ V.map realToFrac evec

getSymmetricEigenpair :: ArpackResults -> Int -> IO (Double, Vector Double)
getSymmetricEigenpair results index = do
  eig <- getSymmetricEigenvalue results index
  evec <- getSymmetricEigenvector results index
  return (eig,evec)

eigsS :: LinearOperator
      -> SWhich
      -> Dimension
      -> HowMany
      -> Tolerance
      -> Iterations
      -> Maybe [(Double,Vector Double)]
eigsS f which n nev tol' iters = unsafePerformIO $ do
  setup <- arpack f (Symmetric which) n nev tol' iters
  results <- parseArpackOutput (Symmetric which) setup
  err <- errCheck results
  if err
    then return Nothing
    else do
      pairs <- mapM (getSymmetricEigenpair results) [0..(nev-1)]
      return (Just pairs)

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

mallocInt :: IO (Ptr CInt)
mallocInt = malloc

ptrInt :: Int -> IO (Ptr CInt)
ptrInt k =
  do ptr <- malloc
     poke ptr (fromIntegral k)
     return ptr

mallocIntArray :: Int -> IO (Ptr CInt)
mallocIntArray = mallocArray

----------

mallocDouble :: IO (Ptr CDouble)
mallocDouble = malloc

ptrDouble :: Double -> IO (Ptr CDouble)
ptrDouble x =
  do ptr <- malloc
     poke ptr (realToFrac x)
     return ptr

mallocDoubleArray :: Int -> IO (Ptr CDouble)
mallocDoubleArray = mallocArray

newVectorDouble :: CInt -> IO (MV.IOVector CDouble)
newVectorDouble = MV.new . fromIntegral

----------

mallocChar :: IO (Ptr CChar)
mallocChar = malloc

mallocCharArray :: Int -> IO (Ptr CChar)
mallocCharArray = mallocArray

ptrCharArray :: String -> IO (Ptr CChar)
ptrCharArray s =
  do ptr <- mallocCharArray (length s + 1)
     pokeArray ptr (map castCharToCChar s)
     return ptr

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

foreign import ccall unsafe "arpack.h dnaupd_"
  c_dnaupd :: Ptr CInt    -- ^ reverse communication flag, to be started with 0
           -> Ptr CChar   -- ^ 'I' for standard problem, 'G' for generalized
           -> Ptr CInt    -- ^ dimension of the eigenproblem
           -> Ptr CChar   -- ^ which eigenpairs to compute
           -> Ptr CInt    -- ^ N : number of eigenpairs to compute
           -> Ptr CDouble -- ^ tolerance
           -> Ptr CDouble -- ^ residual vector
           -> Ptr CInt    -- ^ NCV : number of columns of V matrix
           -> Ptr CDouble -- ^ N x NCV output array
           -> Ptr CInt    -- ^ leading dimension of V
           -> Ptr CInt    -- ^ array for internal parameters, length 11
           -> Ptr CInt    -- ^ pointer of starting locations, length 14
           -> Ptr CDouble -- ^ workd : work array
           -> Ptr CDouble -- ^ workl : work array
           -> Ptr CInt    -- ^ length of workl
           -> Ptr CInt    -- ^ info
           -> IO ()

dnaupd :: ArpackSetup -> IO ()
dnaupd setup =
  c_dnaupd (aIDO setup) (aBMAT setup) (aN setup) (aWHICH setup) (aNEV setup)
           (aTOL setup) (aRESID setup) (aNCV setup) (aV setup) (aLDV setup)
           (aIPARAM setup) (aIPNTR setup) (aWORKD setup) (aWORKL setup)
           (aLWORKL setup) (aINFO setup)

foreign import ccall unsafe "arpack.h dneupd_"
  c_dneupd :: Ptr CInt    -- ^ whether to calculate eigenvectors
           -> Ptr CChar   -- ^ whether to calculate Ritz or Schur vectors
           -> Ptr CInt    -- ^ selects the vectors to be computed
           -> Ptr CDouble -- ^ real part of eigenvalues
           -> Ptr CDouble -- ^ imaginary part of eigenvalues
           -> Ptr CDouble -- ^ approximate eigenvectors
           -> Ptr CInt    -- ^ leading dimension of former array
           -> Ptr CDouble -- ^ real part of the shift
           -> Ptr CDouble -- ^ imaginary part of the shift
           -> Ptr CDouble -- ^ work array
           -> Ptr CChar -> Ptr CInt -> Ptr CChar -> Ptr CInt -> Ptr CDouble
           -> Ptr CDouble -> Ptr CInt -> Ptr CDouble -> Ptr CInt -> Ptr CInt
           -> Ptr CInt -> Ptr CDouble -> Ptr CDouble -> Ptr CInt -> Ptr CInt
           -> IO ()

--------------------------------------------------------------------------------

foreign import ccall unsafe "arpack.h dsaupd_"
  c_dsaupd :: Ptr CInt    -- ^ reverse communication flag, to be started with 0
           -> Ptr CChar   -- ^ 'I' for standard problem, 'G' for generalized
           -> Ptr CInt    -- ^ dimension of the eigenproblem
           -> Ptr CChar   -- ^ which eigenpairs to compute
           -> Ptr CInt    -- ^ N : number of eigenpairs to compute
           -> Ptr CDouble -- ^ tolerance
           -> Ptr CDouble -- ^ residual vector
           -> Ptr CInt    -- ^ NCV : number of columns of V matrix
           -> Ptr CDouble -- ^ N x NCV output array
           -> Ptr CInt    -- ^ leading dimension of V
           -> Ptr CInt    -- ^ array for internal parameters, length 11
           -> Ptr CInt    -- ^ pointer of starting locations, length 14
           -> Ptr CDouble -- ^ workd : work array
           -> Ptr CDouble -- ^ workl : work array
           -> Ptr CInt    -- ^ length of workl
           -> Ptr CInt    -- ^ info
           -> IO ()

dsaupd :: ArpackSetup -> IO ()
dsaupd setup =
  c_dsaupd (aIDO setup) (aBMAT setup) (aN setup) (aWHICH setup) (aNEV setup)
           (aTOL setup) (aRESID setup) (aNCV setup) (aV setup) (aLDV setup)
           (aIPARAM setup) (aIPNTR setup) (aWORKD setup) (aWORKL setup)
           (aLWORKL setup) (aINFO setup)

foreign import ccall unsafe "arpack.h dseupd_"
  c_dseupd :: Ptr CInt    -- ^ whether to calculate eigenvectors
           -> Ptr CChar   -- ^ whether to calculate Ritz or Schur vectors
           -> Ptr CInt    -- ^ selects the vectors to be computed
           -> Ptr CDouble -- ^ eigenvalues
           -> Ptr CDouble -- ^ approximate eigenvectors
           -> Ptr CInt    -- ^ leading dimension of former array
           -> Ptr CDouble -- ^ shift
           -> Ptr CChar -> Ptr CInt -> Ptr CChar -> Ptr CInt -> Ptr CDouble
           -> Ptr CDouble -> Ptr CInt -> Ptr CDouble -> Ptr CInt -> Ptr CInt
           -> Ptr CInt -> Ptr CDouble -> Ptr CDouble -> Ptr CInt -> Ptr CInt
           -> IO ()
