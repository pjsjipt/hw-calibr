### A Pluto.jl notebook ###
# v0.19.6

using Markdown
using InteractiveUtils

# ╔═╡ 53a2a9e7-760a-48b3-877d-dac8e878a6da
begin
	import Pkg
	Pkg.activate(".")
	using Polynomials
	using GLMakie
	using ConstructiveGeometry
end
	

# ╔═╡ aac6b8d1-11fb-4f99-ae99-a4399d8f9c07
using FileIO


# ╔═╡ 0de5472e-e26c-11ec-1beb-15a78a5cb13b
md"""
# Contração de túnel de vento

Vamos usar uma cúbica
"""

# ╔═╡ 00c86c4b-8e46-4d83-8d6e-3a77ed3f99ae


# ╔═╡ 145112a4-5a96-48c3-8d17-a79edff2bb5f
function contrac_profile(L, H, b, ang=2)
	α = -ang * π / 180   
	H₁ = H
	H₂ = b*H
	ΔH = (H₂ - H₁) / 2

	a₀ = -H₂ / 2  # Center line is the symmetry axis
	a₁ = 0.0
	a₂ = (3ΔH/L - α) / L
	a₃ = (α - 2ΔH/L) / L^2

	return Polynomial([a₀, a₁, a₂, a₃])
end


# ╔═╡ 41167b23-5cd6-41cd-8f83-fc9fee273e41
begin
	H = 20 # Altura da seção de testes
	b = 3
	L = 3*H # Comprimento da contração
	H2 = b*H
	esp = 3.0 # Espessura da parede
	A = 12.0 # Largura das flanges
	p = contrac_profile(L, H, b, 2)
	p1 = p - Polynomial([esp])
	d = 4.0 # Diâmetro dos furos
	Dm =25.4 * 3/4
	Lm = 1.5*L
	
end

# ╔═╡ 7e589de0-aaa3-4db0-ac61-21d5698e4381
function make_flange(Hi, nh, esp, A=12.0, dh=4.0, lh=6.0)
	Le = Hi + 2*(esp+A)
	Li = Hi + esp

	c1 = translate([-Le/2, -Le/2, 0])*cube(Le, Le, esp)
	c2 = translate([-Li/2, -Li/2, -esp/2]) * cube(Li,Li,2*esp)
	cyl = translate([0, 0, -2*esp]) * cylinder(6*esp, d/2)
	Xf = Le/2-lh   # Distance of the holes from flange center
	# Start with external cube minus internal cube and the corner holes
	flange = c1 - c2 - translate([-Xf, -Xf, 0]) * cyl - translate([Xf, -Xf, 0]) * cyl - translate([Xf, Xf, 0]) * cyl - translate([-Xf, Xf, 0]) * cyl

	# If n > 2 get the intermediate wholes
	if nh > 2
		Xi = range(-Xf, Xf, length=nh)
		for i in 2:nh-1
			flange = flange - translate([Xi[i], -Xf, 0]) * cyl - translate([Xi[i], +Xf, 0]) * cyl - translate([-Xf, Xi[i], 0]) * cyl - translate([+Xf, Xi[i], 0]) * cyl
		end
	end

	return flange
end


# ╔═╡ aea743f9-966e-4e9b-9136-0c4cf28ffbb9
function make_contrac(np, H1, H2, L, esp, ang)
	xc = range(0, stop=L, length=np)
	δ = xc[2]
	p = contrac_profile(L, H1, H2/H1, ang)
	p1 = p - Polynomial([esp])
	
	y0 = p.(xc)
	y1 = p1.(xc)

	xx0 = [-δ; xc; L+δ; L+δ; reverse(xc); -δ]
	yy0 = [y0[1]; y0; y0[end]; -y0[end]; reverse(-y0); -y0[1]]
	xx1 = [xc; reverse(xc)]
	yy1 = [y1; reverse(-y1)]
	
	poly0 = polygon([pt for pt in zip(xx0, yy0)])
	poly1 = polygon([pt for pt in zip(xx1, yy1)])
	solid0 = rotate((0,-90,0))*translate([0,0,-H2/2])*linear_extrude(H2) * poly0
	solid1 = rotate((0,-90,0))*translate([0,0,-H2/2-esp])*linear_extrude(H2+2*esp) * 
				poly1	
	contrac_int = solid0 ∩ (rotate(90)*solid0)
	contrac_ext = solid1 ∩ (rotate(90)*solid1)	

	return contrac_ext - contrac_int
end


# ╔═╡ edd6b1aa-df65-4c56-b8d5-9269ba7228de


# ╔═╡ 1ce1a48b-a2e3-4a32-a153-1a7d915fc094
plot(make_contrac(100, 20, 100, 100, 3, 0))

# ╔═╡ b469ac94-4043-4edc-a7ed-8f5de266589f
plot(make_flange(20, 3, esp))

# ╔═╡ a9f14606-9ebf-49cd-9ba3-6bf8c6026bca
let
	x = range(0.0, stop=L, length=100)
	y = p.(x)
	y1 = p1.(x)
	fig = Figure()
	ax = Axis(fig[1,1], aspect=DataAspect())
	lines!(ax, x, y)
	lines!(ax, x, -y)
	lines!(ax, x, y1)
	lines!(ax, x, -y1)
	fig
end

# ╔═╡ d4bd380f-a762-45ff-8e11-006f8610c42f
md"""
## Criar a geometria 3d 

"""

# ╔═╡ 96700874-1b3d-4e24-a8b1-607b69fad6ed
begin
	np = 100
	xc = range(0, stop=L, length=np)
	δ = xc[2]
	y0 = p.(xc)
	y1 = p1.(xc)

	xx0 = [-δ; xc; L+δ; L+δ; reverse(xc); -δ]
	yy0 = [y0[1]; y0; y0[end]; -y0[end]; reverse(-y0); -y0[1]]
	xx1 = [xc; reverse(xc)]
	yy1 = [y1; reverse(-y1)]
	
	poly0 = polygon([pt for pt in zip(xx0, yy0)])
	poly1 = polygon([pt for pt in zip(xx1, yy1)])
	nothing
end


# ╔═╡ 296a9294-2e68-4f4c-96ee-af3de99b74d6
solid0 = rotate((0,-90,0))*translate([0,0,-b*H/2])*linear_extrude(b*H) * poly0

# ╔═╡ bcc6c25f-ce56-4ca0-8ba1-9963672e94e8
solid1 = rotate((0,-90,0))*translate([0,0,-b*H/2-esp])*linear_extrude(b*H+2*esp) * poly1

# ╔═╡ 85d7c4b7-40de-4f44-9ad8-1caa7dfaae8f
plot(solid1)

# ╔═╡ f8999f20-a3fc-49a2-8aef-7d88415911ca
contrac_int = solid0 ∩ (rotate(90)*solid0)

# ╔═╡ 39b5e912-42a2-49c4-897d-1d67e67d41e1
contrac_ext = solid1 ∩ (rotate(90)*solid1)

# ╔═╡ 5d062dc0-b35f-42fd-b463-f0d5d445889d
contrac = setdiff(contrac_ext, contrac_int)

# ╔═╡ 4840bd2e-ac34-44e6-a0c5-54dac7996b75
begin
	plot(contrac)
end

# ╔═╡ d56286f2-66b4-4ca9-8e83-7ede993f97b7
# Aba da seção de saída
aba0 = let
	Le = H + 2*(esp+A)
	Li = H + esp

	c1 = translate([-Le/2, -Le/2, 0])*cube(Le, Le, esp)
	c2 = translate([-Li/2, -Li/2, -esp/2]) * cube(Li,Li,2*esp)
	cyl = translate([0, 0, -2*esp]) * cylinder(6*esp, d/2)
	Xf = Le/2-1.5*d
	
	c1-c2 - translate([-Xf, -Xf, 0]) * cyl - translate([0, -Xf, 0]) * cyl - translate([+Xf, -Xf, 0]) * cyl - translate([-Xf, 0, 0]) * cyl - translate([-Xf, +Xf, 0]) * cyl - translate([0, +Xf, 0]) * cyl - translate([+Xf, +Xf, 0]) * cyl - translate([+Xf, 0, 0]) * cyl 
end

# ╔═╡ 2ac01d6a-d91e-46d3-99ab-48222144e59c
plot(aba0)

# ╔═╡ e68c6bdb-3323-4549-b7a2-e958695d33f8
# Aba externa
aba1 = let
	Le = H2 + 2*(esp+A)
	Li = H2 + esp

	c1 = translate([-Le/2, -Le/2, 0])*cube(Le, Le, esp)
	c2 = translate([-Li/2, -Li/2, -esp/2]) * cube(Li,Li,2*esp)
	cyl = translate([0, 0, -2*esp]) * cylinder(6*esp, d/2)
	Xf = Le/2-1.5*d
	
	c1-c2 - translate([-Xf, -Xf, 0]) * cyl - translate([0, -Xf, 0]) * cyl - translate([+Xf, -Xf, 0]) * cyl - translate([-Xf, 0, 0]) * cyl - translate([-Xf, +Xf, 0]) * cyl - translate([0, +Xf, 0]) * cyl - translate([+Xf, +Xf, 0]) * cyl - translate([+Xf, 0, 0]) * cyl - translate([+Xf/2, +Xf, 0]) * cyl - translate([-Xf/2, +Xf, 0]) * cyl - translate([+Xf/2, -Xf, 0]) * cyl - translate([-Xf/2, -Xf, 0]) * cyl - translate([-Xf, -Xf/2, 0]) * cyl - translate([-Xf, +Xf/2, 0]) * cyl - translate([+Xf, -Xf/2, 0]) * cyl - translate([+Xf, +Xf/2, 0]) * cyl 
end	

# ╔═╡ bba1b9dc-185b-497f-ac8b-0b34f2a27985
wtcontrac = contrac ∪ (translate([0,0,L-esp])*aba0) ∪ aba1

# ╔═╡ 6f00d252-af67-4817-bdb6-dfbad89ec03d
plot(wtcontrac)

# ╔═╡ 730868ab-3f1a-44da-bf99-52e122116031
caixa = let
	Li = H2
	Le = H2+2*esp
	W = H2
	c1 = translate([-Li/2, -Li/2, -W/2]) * cube(Li,Li,2*W)
	c2 = translate([-Le/2, -Le/2, 0]) * cube(Le,Le,W)
	(c2-c1) ∪ aba1 ∪ (translate([0,0,W-esp])*aba1)
end

# ╔═╡ c50d0bae-7d62-4dd0-9596-9c82fdb112a1
plot(caixa)

# ╔═╡ a6d4b902-b7df-4f5a-aa5a-de5426e790c3
caixa2 = let
	Li = H2
	Le = H2+2*esp
	W = H2 / 2
	c1 = translate([-Li/2, -Li/2, -W/2]) * cube(Li,Li,2*W)
	c2 = translate([-Le/2, -Le/2, 0]) * cube(Le,Le,W)
	(c2-c1) ∪ aba1 ∪ (translate([0,0,W-esp])*aba1)
end

# ╔═╡ 994246b5-78c6-4909-942d-57521d2e2e8e
plot(caixa2)

# ╔═╡ af66294c-a432-4d71-94b3-cf830a7e611a
md"""
## Expansão de entrada
"""

# ╔═╡ 5b3434df-d909-4598-b4b4-47dfc880ba1b
entrada = let
	# make_contrac(np, H1, H2, L, esp, ang)
	Hm = Dm + 2*esp
	ee = make_contrac(100, Dm, H2, Lm, esp, 0)
	cf = translate([-Hm/2, -Hm/2, Lm-1]) * cube(Hm, Hm, Hm)
	furo = translate([0,0,Lm-Hm]) * cylinder(3*Hm, Dm/2)
	aba1 ∪ ee ∪ (cf - furo)
		
	#=
	Hm = Dm + 2*esp
	Le = H2 + 2*esp
	Li = H2
	Lc = Le*Lm / (Le - Hm)
	si = translate([-Li/2, -Li/2]) * square(Li,Li)
	se = translate([-Le/2, -Le/2]) * square(Le,Le)
	s = se - si
	ce = translate([-Le/2, -Le/2, Lm]) * cube(Le,Le,Lc)
	redu = cone([0,0,Lc]) * s - ce
	cf = translate([-Hm/2, -Hm/2, Lm-1]) * cube(Hm, Hm, Hm)
	furo = translate([0,0,Lm-Hm]) * cylinder(3*Hm, Dm/2)
	aba1 ∪ redu ∪ (cf - furo)
	=#
end

# ╔═╡ 66b347c9-fa24-4411-8e6b-0078e58bb94b
let

end


# ╔═╡ 22fe12a7-1532-41e6-b23b-a6023b8c627b
plot(entrada)

# ╔═╡ ce41e94a-318c-4300-8200-f875e4773585
md"""
## Base para colocar o anemômetro
"""

# ╔═╡ ba7f8d3c-7116-499f-8afd-a6601ed63888
base = let
	dd = 10
	rr = dd/2
	Wb = H + 2*esp 
	base  = translate([-Wb/2, -esp-H/2, 0]) * cube(Wb, esp, H+esp)
	cyl = 	translate([0, -H/2-2*esp, H/2+esp]) * (rotate((0,0,-90)) * cylinder(5*esp, rr))
	cc = translate([-dd/2, -H/2 - 2*esp, esp+H/2])*cube(dd, 5*esp, H)
	aba0 ∪ (base - cyl - cc)
end


# ╔═╡ 099cff15-e1ac-4150-8ce0-9beb0b631095
plot(base)

# ╔═╡ 7dcaad52-cfc0-4f1b-aa7f-303276077363
md"""

## Vamos salvar tudo

"""

# ╔═╡ 0824743a-6c1b-4e43-b690-49ee70a034af
begin
	save("wtcontrac.stl", wtcontrac)
	save("caixa.stl", caixa)
	save("caixa2.stl", caixa2)
	save("entrada.stl", entrada)
	save("base.stl", base)
end

# ╔═╡ 3d50dd45-2d31-4070-b819-e4864d279b53
typeof(base)

# ╔═╡ Cell order:
# ╠═0de5472e-e26c-11ec-1beb-15a78a5cb13b
# ╠═53a2a9e7-760a-48b3-877d-dac8e878a6da
# ╠═00c86c4b-8e46-4d83-8d6e-3a77ed3f99ae
# ╠═145112a4-5a96-48c3-8d17-a79edff2bb5f
# ╠═41167b23-5cd6-41cd-8f83-fc9fee273e41
# ╠═7e589de0-aaa3-4db0-ac61-21d5698e4381
# ╠═aea743f9-966e-4e9b-9136-0c4cf28ffbb9
# ╠═edd6b1aa-df65-4c56-b8d5-9269ba7228de
# ╠═1ce1a48b-a2e3-4a32-a153-1a7d915fc094
# ╠═b469ac94-4043-4edc-a7ed-8f5de266589f
# ╠═2ac01d6a-d91e-46d3-99ab-48222144e59c
# ╠═a9f14606-9ebf-49cd-9ba3-6bf8c6026bca
# ╠═d4bd380f-a762-45ff-8e11-006f8610c42f
# ╠═96700874-1b3d-4e24-a8b1-607b69fad6ed
# ╠═296a9294-2e68-4f4c-96ee-af3de99b74d6
# ╠═bcc6c25f-ce56-4ca0-8ba1-9963672e94e8
# ╠═85d7c4b7-40de-4f44-9ad8-1caa7dfaae8f
# ╠═f8999f20-a3fc-49a2-8aef-7d88415911ca
# ╠═39b5e912-42a2-49c4-897d-1d67e67d41e1
# ╠═5d062dc0-b35f-42fd-b463-f0d5d445889d
# ╠═4840bd2e-ac34-44e6-a0c5-54dac7996b75
# ╠═d56286f2-66b4-4ca9-8e83-7ede993f97b7
# ╠═e68c6bdb-3323-4549-b7a2-e958695d33f8
# ╠═bba1b9dc-185b-497f-ac8b-0b34f2a27985
# ╠═6f00d252-af67-4817-bdb6-dfbad89ec03d
# ╠═730868ab-3f1a-44da-bf99-52e122116031
# ╠═c50d0bae-7d62-4dd0-9596-9c82fdb112a1
# ╠═a6d4b902-b7df-4f5a-aa5a-de5426e790c3
# ╠═994246b5-78c6-4909-942d-57521d2e2e8e
# ╠═af66294c-a432-4d71-94b3-cf830a7e611a
# ╠═5b3434df-d909-4598-b4b4-47dfc880ba1b
# ╠═66b347c9-fa24-4411-8e6b-0078e58bb94b
# ╠═22fe12a7-1532-41e6-b23b-a6023b8c627b
# ╠═ce41e94a-318c-4300-8200-f875e4773585
# ╠═ba7f8d3c-7116-499f-8afd-a6601ed63888
# ╠═099cff15-e1ac-4150-8ce0-9beb0b631095
# ╠═7dcaad52-cfc0-4f1b-aa7f-303276077363
# ╠═aac6b8d1-11fb-4f99-ae99-a4399d8f9c07
# ╠═0824743a-6c1b-4e43-b690-49ee70a034af
# ╠═3d50dd45-2d31-4070-b819-e4864d279b53
