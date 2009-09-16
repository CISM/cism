EXTRA_DIST = index.html template.config error.cpt

test_b = test_B_non-lin_50km_10.00a.nc test_B_lin_50km_10.00a.nc \
         test_B_non-lin_20km_2.00a.nc  test_B_lin_20km_2.00a.nc \
         test_B_non-lin_10km_0.50a.nc  test_B_lin_10km_0.50a.nc \
         test_B_non-lin_5km_0.10a.nc   test_B_lin_5km_0.10a.nc 
test_c = test_C_non-lin_50km_5.00a.nc test_C_lin_50km_5.00a.nc \
         test_C_non-lin_20km_1.00a.nc test_C_lin_20km_1.00a.nc \
         test_C_non-lin_10km_0.50a.nc test_C_lin_10km_0.50a.nc \
         test_C_non-lin_5km_0.10a.nc  test_C_lin_5km_0.10a.nc 
test_d = test_D_non-lin_50km_5.00a.nc test_D_lin_50km_5.00a.nc \
         test_D_non-lin_20km_1.00a.nc test_D_lin_20km_1.00a.nc

all_tests = $(test_b) $(test_c) $(test_d)

plots:		vole_c.ps vole_d.ps errors.ps time.ps
png:		vole_c.png vole_d.png errors.png time.png $(subst nc,png,$(all_tests))

errors.ps:	$(all_tests)
		$(PYTHON) $(top_srcdir)/src/python/plot_errors.py $^ $@

time.ps:	data
		$(PYTHON) $(top_srcdir)/src/python/plot_time.py $^ $@

vole_b.ps:	$(test_b)
		$(PYTHON) $(top_srcdir)/src/python/plot_vol.py $^ $@

vole_c.ps:	$(test_c)
		$(PYTHON) $(top_srcdir)/src/python/plot_vol.py $^ $@

vole_d.ps:	$(test_d)
		$(PYTHON) $(top_srcdir)/src/python/plot_vol.py $^ $@

data:		$(all_tests)

%.ps:	%.nc
	$(PYTHON) $(top_srcdir)/src/python/plot_diff.py -T-1 $< $@

%.nc:	$(top_srcdir)/tests/exact-isothermal/tests/template.config
	$(PYTHON) $(top_srcdir)/tests/exact-isothermal/scripts/run_verif.py -f $* $<

%.png:	%.ps
	$(CONVERT) -trim $< $@
