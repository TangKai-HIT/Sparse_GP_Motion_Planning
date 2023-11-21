function constraints = gpmp_cons_init(posRange, velRange, accRange, jerkRange)
%GPMP_CONS_INIT Get a default struct of gpmp constraints
%   posRange, velRange, accRange, jerkRange: dim X 2, e.g. [[min_x, max_x]; [min_y, max_y]; ... ]

constraints.posRange = posRange;
constraints.velRange = velRange;
constraints.accRange = accRange;
constraints.jerkRange = jerkRange;
end

