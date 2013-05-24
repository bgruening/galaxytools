# encoding: utf-8
#
# Jekyll category pagination.
# https://github.com/MrWerewolf/jekyll-category-pagination
#
# Copyright (c) 2012 Ryan Seto <mr.werewolf@gmail.com>
# Licensed under the MIT license (http://www.opensource.org/licenses/mit-license.php)
#
# Place this script into the _plugins directory of your Jekyll site.
#
module Jekyll

  # Change the Page object so it puts the path in the url.  This makes the dir
  # method work as intended.
  class Page
    # The template of the permalink.
    #
    # Returns the template String.
    def template
      if self.site.permalink_style == :pretty && !index? && html?
        "/:path/:basename/"
      else
        "/:path/:basename:output_ext"
      end
    end

    # The generated relative url of this page. e.g. /about.html.
    #
    # Returns the String url.
    def url
      return @url if @url

      url = if permalink
        permalink
      else
        {
          "path"       => @dir,
          "basename"   => self.basename,
          "output_ext" => self.output_ext,
        }.inject(template) { |result, token|
          result.gsub(/:#{token.first}/, token.last)
        }.gsub(/\/\//, "/")
      end

      # sanitize url
      @url = url.split('/').reject{ |part| part =~ /^\.+$/ }.join('/')
      @url += "/" if url =~ /\/$/
      @url
    end

    # Convert this Page's data to a Hash suitable for use by Liquid.
    #
    # Returns the Hash representation of this Page.
    def to_liquid
      self.data.deep_merge({
        "url"        => self.url,
        "content"    => self.content })
    end

    # Obtain destination path.
    #
    # dest - The String path to the destination dir.
    #
    # Returns the destination file path String.
    def destination(dest)
      # The url needs to be unescaped in order to preserve the correct
      # filename.
      path = File.join(dest, CGI.unescape(self.url))
      path = File.join(path, "index.html") if self.url =~ /\/$/
      path
    end
  end

  class Pagination < Generator
    # This generator is safe from arbitrary code execution.
    safe true

    # Paginates the blog's posts. Renders the index.html file into paginated
    # directories, e.g.: page2/index.html, page3/index.html, etc and adds more
    # site-wide data.
    #
    # site - The Site.
    # page - The index.html Page that requires pagination.
    #
    # {"paginator" => { "page" => <Number>,
    #                   "per_page" => <Number>,
    #                   "posts" => [<Post>],
    #                   "total_posts" => <Number>,
    #                   "total_pages" => <Number>,
    #                   "previous_page" => <Number>,
    #                   "next_page" => <Number> }}
    def paginate(site, page)
      category = page.dir.split('/').last
      if category != nil and site.site_payload['site']['categories'].has_key? category
        # If we're in a category's folder, paginate by category.
        all_posts = site.site_payload['site']['categories'][category]
      else
        all_posts = site.site_payload['site']['posts']
      end

      pages = Pager.calculate_pages(all_posts, site.config['paginate'].to_i)
      (1..pages).each do |num_page|
        pager = Pager.new(site.config, num_page, all_posts, pages)
        if num_page > 1
          newpage = Page.new(site, site.source, page.dir, page.name)
          newpage.pager = pager
          newpage.dir = File.join(page.dir, "page#{num_page}")
          site.pages << newpage
        else
          page.pager = pager
        end
      end
    end
  end
end
