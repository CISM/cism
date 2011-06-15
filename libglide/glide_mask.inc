#ifndef GLIDE_MASK_INC

!Bits 1:0 - Ice presence (0 if no ice)
!NOTE: If bit 2 is activated, bit 1 must be activated
#define GLIDE_MASK_HAS_ICE             1
#define GLIDE_MASK_THIN_ICE            3
#define GLIDE_ICE_PRESENCE_BITS        3

!Bits 4:2 - Type of base (Land or ocean - grounding line has both bits on)
!The 16 bit specifies grounding line - I am setting this up so that those points are treated as grounded
!as well.
#define GLIDE_MASK_LAND                4
#define GLIDE_MASK_OCEAN               8
#define GLIDE_MASK_GROUNDING_LINE      16 
#define GLIDE_BASE_TYPE_BITS           28

!Bit 5: Identifies a margin (jump from zero to nonzero thickness)
!The type of margin is determined by whether the ice here is grounded or floating.
#define GLIDE_MASK_MARGIN 32

!Bit 6: Identifies a dirichlet condition
!The velocity at points marked with this bit should be held constant.
#define GLIDE_MASK_DIRICHLET_BC 64

!Bit 7: Identifies a computational domain boundary.  These are normally just activated on the edges of the domain,
!unless there is a domain decomposition (in which case they may be missing)
#define GLIDE_MASK_COMP_DOMAIN_BND 128 

!Checks for an iceless square

!Checks for open ocean with no ice
#define GLIDE_IS_OCEAN(mask)       (iand(mask, GLIDE_MASK_OCEAN) == GLIDE_MASK_OCEAN .and. GLIDE_NO_ICE(mask))

!Checks for land with no ice
#define GLIDE_IS_LAND(mask)        (iand(mask, GLIDE_MASK_LAND) == GLIDE_MASK_LAND .and. GLIDE_NO_ICE(mask))

!Checks for the presence of any ice, dynamic or not
#define GLIDE_HAS_ICE(mask)        (iand(mask, GLIDE_MASK_HAS_ICE) == GLIDE_MASK_HAS_ICE)

!Checks for a lack of ice
#define GLIDE_NO_ICE(mask)         (iand(mask, GLIDE_MASK_HAS_ICE) == 0)

!Checks for the presence of ice that is below the ice dynamics limit
#define GLIDE_IS_THIN(mask)          (iand(mask,GLIDE_MASK_THIN_ICE)  == GLIDE_MASK_THIN_ICE)

!Checks for any ice, dynamic or not, that is on an ice shelf.
#define GLIDE_IS_FLOAT(mask)         (iand(mask,GLIDE_MASK_OCEAN) == GLIDE_MASK_OCEAN .and. GLIDE_HAS_ICE(mask))

!Checks for any ice, dynamic or not, that is grounded
#define GLIDE_IS_GROUND(mask)      (iand(mask,GLIDE_MASK_LAND) == GLIDE_MASK_LAND .and. GLIDE_HAS_ICE(mask))

!Checks for any ice, dynamic or not, that is on the grounding line
#define GLIDE_IS_GROUNDING_LINE(mask) (iand(mask, 17) == 17)

!Checks for any ice, dynamic or not, that is either floating *or* on the grounding line
#define GLIDE_IS_FLOAT_OR_GNDLINE(mask) (iand(mask, 24) > 0 .and. GLIDE_HAS_ICE(mask))

!Checks for any ice, dynamic or not, that is either grounded *or* on the grounding line
#define GLIDE_IS_GROUND_OR_GNDLINE(mask) (iand(mask, 20) > 0 .and. GLIDE_HAS_ICE(mask))

!Checks whether this is an ice margin (thickness jumps from 0 to non-zero at this point)
#define GLIDE_IS_MARGIN(mask) (iand(mask, GLIDE_MASK_MARGIN) == GLIDE_MASK_MARGIN)

!Checks whether this is a margin in contact with the ocean, floating or not
#define GLIDE_IS_MARINE_ICE_EDGE(mask)    (GLIDE_IS_MARGIN(mask) .and. GLIDE_IS_FLOAT_OR_GNDLINE(mask))

!Checks whether this is a margin in contact with the ocean
!41 = 32+8+1
#define GLIDE_IS_CALVING(mask)            (iand(mask, 41) == 41)

!Checks whether this is a land margin
!37 = 32 + 4 + 1
#define GLIDE_IS_LAND_MARGIN(mask)        (iand(mask, 37) == 37)

!Checks whether a dirichlet boundary has been defined at this point
#define GLIDE_IS_DIRICHLET_BOUNDARY(mask) (iand(mask, GLIDE_MASK_DIRICHLET_BC) == GLIDE_MASK_DIRICHLET_BC)

!Checks whether we are at the edge of the computational domain *and* there is ice in this square
#define GLIDE_IS_COMP_DOMAIN_BND(mask)    (iand(mask, 129) == 129)


#endif

